//#################################******************************##########################/////////
//#################################***Simulation Paramaeters***##########################/////////
//#################################******************************##########################/////////
//Reynolds Number - 0.4
//Contact angle -60
//Interface Refinement-10
//Global Refinement- 10-3=7  
//Phase A - air 
//Phase B - Water 
//Updated - Yes
//#################################******************************##########################/////////
//#################################******************************##########################/////////
//This is a simualtion to visualize the flow field near a moving contact line. 
//The geometry of the problem is a sqaure domain of size 15 x 15 (ALL LENGTHS IN mm )
//it is a 15 x 15 square with interface in the middle, horizontally. 
//Author- Anvesh 
//The centre of the domain is at the centre of the left wall. 
//We are working in SI units. 
//Date -13-April-2022
//Comments: 
//Status : working 
/* FOR COMPILING :
  qcc   -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
*/
//Libraries used - 

//#include "navier-stokes/conserving.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "contact.h"
#include "tension.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "output_vtu_foreach.h"

int maxlevel = 9;              // Maximum mesh refinement
char name_vtk[100];             // vtk file name decleration.
double U0;
double H0;


        #define grav  -9.81 // gravitational acceleration
        #define rhoL 996 //density of water
        #define muL 0.00089 //viscosity of water
        #define surf  0.072  // surface tension air-water
        #define rhoG 1.2 //density of air
        #define muG  0.0000181 // viscosity of air
        #define lc 2.7e-3// capillary length 
	#define T_end 100 
        #define Uplate  -0.00013145 // plate velocity 
        #define f_tol 0.1    // The tolerance given to the vof field f. 
        #define ux_tol 0.05    // The tolerance given to the ux 
        #define uy_tol 0.05    // The tolerance given to the uy. 
//From here onwards we define the 9 constants for the 8 degree polynomial we are
//fitting for the initial meniscus shape from the 
//final steady state shape from earlier simulations


double h0;

vector h[];  //HEIGHT FUNCTION 
double theta0;

//make sure that the boundary conditions for the face-centered velocity field are consistent with the centered velocity field (this affects the advection term).
uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

int main()
{
        L0 = 0.015;   // Size of the square box -- Upon checking where the interface becomes flat
	origin (0, -L0/2);  // Origin is at the bottom centre of the box
	N = 1024;
        stokes = true;
        f.sigma = surf;
        f.height = h;
        display_control (maxlevel, 6, 15);

        theta0 = 120*pi/180.0;
        h.t[left] = contact_angle (theta0); // Left contact angle near the moving wall 
        h.t[right] = contact_angle (pi/2);  // right contact angle of 90 degrees. 

        //The viscosity and desinties of the two fluids is specified here. 
        rho2 = rhoL;   // fluid 2 is given by f = 0.
        mu2 = muL;
        rho1 = rhoG;   // fluid 1 is given by f =1. 
        mu1 = muG;        
	
	
	run();

}


event init (t = 0)
{
	if(!restore (file = "---- "))
		{
	//Here the approximate static meniscus shape is given as an initial condition.  
//the top fluid has f = 0 and is gas and the bottom fluid is f =1 and is liquid. 
//refer: http://basilisk.fr/src/two-phase.h

        fraction (f,  0.001557050 + y + 0.0027/(tan(theta0)*exp((x)/0.0027)));
	boundary ({f});
/*
		foreach()
		{
		
			foreach_dimension()
			{
			u.x[]=0.0;

			}
		
		}
*/
		}

}



// gravity is given in the vertically down direction.(-9.81)
event acceleration (i++
{
        face vector av = a;
        foreach_face(x)
                av.y[] = grav;
}

// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(Uplate);


u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.0);



// Printing out standard text outputs on the screen
event logfile (i+=50)
        fprintf (stderr, "%d %g\n", i, t);


// The interface profile is extracted for convergence check 
// the reference code is taken from: http://basilisk.fr/Miguel/spreading.c 
event profile(t+=0.1   ; t <= T_end ) 
{
  char int_prof[80];
  sprintf(int_prof,"interface_profile_t%2f.dat", t);
  FILE * fp1 = fopen(int_prof,"w");
  output_facets (f, fp1);
}

/*
void backup_fields(scalar f,vector u,int nf)
{
    char name[80], subname[80];
    FILE *fp;
    nf > 0 ? sprintf(name,"sol-%4.4d_n%3.3d.vtu",nf,pid()) : sprintf(name,"sol-0_n%3.3d.vtu",pid());
    fp = fopen(name,"w");
    output_vtu_ascii_foreach((scalar *){f},(vector *){u},N,fp,false);
    fclose(fp);

    #if _MPI
        if (pid() == 0)
        {
            nf > 0 ? sprintf(name,"sol-%d.pvtu",nf) : sprintf(name,"sol-0.pvtu");
            nf > 0 ? sprintf(subname,"sol-%4.4d",nf) : sprintf(subname,"sol-0");
            fp = fopen(name,"w");
            output_pvtu_ascii((scalar *){f},(vector *){u},N,fp,subname);
            fclose (fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
}

*/
event dumpfile(t = 0; t+= 0.01; t <= T_end)
{
    char name[80];
    sprintf(name,"dump-%g",t*100);
    scalar pid[];
    foreach()
    {
        pid[] = fmod(pid()*(npe() + 37),npe());
    }
    boundary({pid});
    dump(name);
}


		
event adapt (i++) {
  adapt_wavelet ((scalar*){f,u}, (double[]){f_tol,ux_tol,uy_tol},maxlevel , maxlevel-3 );
}
/*
event logfile (t = 0; t += 0.01; t <= T_end)
{
    if (pid() == 0)
    {
        FILE *trackfile;
        trackfile = fopen("track.out","aw");
        fprintf (trackfile,"%.4e\t %.4e\n",t,dt);
        fclose(trackfile);
    }  
 
    static int nf = 0;
    backup_fields(f,u,nf);
    nf++;
}
*/
