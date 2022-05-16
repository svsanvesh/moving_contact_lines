//This is a simualtion to visualize the flow field near a moving contact line. 
//The geometry of the problem is a sqaure domain of size L=5*l_c ; where l_c = 3. (ALL LENGTHS IN m )
//it is a 15x15 square with interface in the middle, horizontally. 
//Author- Anvesh 
//The centre of the domain is at the centre of the left wall. 
//We are working in SI units. 
//Date - 2-Feb-2022
//
//Comments: 
//Status : working 
/* FOR COMPILING :
  qcc   -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
*/
//
//Libraries used - 

//#include "navier-stokes/conserving.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "adapt_wavelet_leave_interface.h"
#include "contact.h"
#include "tension.h"
#include "two-phase.h"
#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f)  (clamp(f,0,1)*(mu1 - mu2) + mu2)

int maxlevel = 9;              // Maximum mesh refinement
char name_vtk[100];             // vtk file name decleration.
double U0;
double H0;


        #define grav  -9.81 // gravitational acceleration
        #define rhoL 940  //density of water
        #define muL 0.0088 //viscosity of water
        #define surf  0.021  // surface tension air-water
        #define rhoG 1.2 //density of air
        #define muG  0.0000181 // viscosity of air
        #define lc 1.51e-3// capillary length 
	#define T_end 10
        #define Uplate  -0.0001 // plate velocity 
        #define f_tol 0.0001    // The tolerance given to the vof field f.
        #define ux_tol 0.05    // The tolerance given to the ux
        #define uy_tol 0.05    // The tolerance given to the uy.

//From here onwards we define the 9 constants for the 8 degree polynomial we are
//fitting for the initial meniscus shape from the 
//final steady state shape from earlier simulations

	#define p1  9.3395e+13 
	#define p2  -5.3835e+11
	#define p3 -6.1069e+09
	#define p4 2.0621e+07
	#define p5 2.8029e+05 
	#define p6 -1.6932e+03
	#define p7 8.148
	#define p8 -0.0486
	#define p9 -0.0004953

double h0;

vector h[];  //HEIGHT FUNCTION 
double theta0 ;

//make sure that the boundary conditions for the face-centered velocity field are consistent with the centered velocity field (this affects the advection term).
uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

int padding=6;
int main()
{
        L0 = 0.015;            // Size of the square box
        U0 = -0.001 ;             // Velocity of the left plate
	origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
	N = 128;
//        stokes = true;
        f.sigma = surf;
        f.height = h;
        display_control (maxlevel, 6, 15);

        theta0 = 120*pi/180.0;
        h.t[top] = contact_angle (theta0); // Left contact angle near the moving wall 
        h.t[bottom] = contact_angle (pi/2);  // right contact angle of 90 degrees. 

        //The viscosity and desinties of the two fluids is specified here. 
        rho2 = rhoL;   // fluid 2 is given by f = 0.
        mu2 = muL;
	rho1 = rhoG;   // fluid 1 is given by f =1. 
	mu1 = muG;
        run();

}


event init (t = 0)
{
//Here the approximate static meniscus shape is given as an initial condition.  
//the top fluid has f = 0 and is gas and the bottom fluid is f =1 and is liquid. 
//refer: http://basilisk.fr/src/two-phase.h

//        fraction (f,y-( p1*x*x*x*x*x*x*x*x + p2*x*x*x*x*x*x*x + p3*x*x*x*x*x*x + p4*x*x*x*x*x + p5*x*x*x*x + p6*x*x*x + p7*x*x + p8*x + p9));
//        f.refine = f.prolongation = fraction_refine;
        fraction (f,  ( x + lc/(tan(theta0)*exp((-y+ 0.0075)/lc))));
        boundary ({f});
}



// gravity is given in the vertically down direction.(-9.81)
event acceleration (i++)
{
        face vector av = a;
        foreach_face(x)
                av.x[] = grav;
}

// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.0);


u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);


u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(Uplate);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.0);



// Printing out standard text outputs on the screen
event logfile (i+=50)
        fprintf (stderr, "%d %g\n", i, t);


// The interface profile is extracted for convergence check 
// the reference code is taken from: http://basilisk.fr/Miguel/spreading.c 


// Both dump files and bviewable files are generated here. 
event dumpfile (t = 0.;t += 0.01;t <= T_end)
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


// files generated for viewing in bview

        sprintf (name, "bview-%g",t*100);
        dump (name);

}


/*
event videos ( t+=0.001   ; t <= T_end )
{

        output_ppm (f, file = "f_plate_adv.mp4",8192,
                        min = 0, max = 1.0, linear = true);
//      This snippet of code help put time on top right corner. 
//      reference: http://basilisk.fr/src/examples/breaking.c
        char fname[100];
        sprintf (fname, " t = %.6f ", t );
        draw_string (fname, pos=2, size = 60);
        squares("f",min = 0, max = 1.0, linear = true);
        cells();
        draw_vof ("f" );
        save ("fd.mp4");
        clear();
}
*/
/*
//Here the code makes sure the refinement of the interface is high. 
event adapt (i += 5) {
  adapt_wavelet ((scalar*){f}, (double[]){f_tol},maxlevel);
}
*/

event adapt (i += 5) {
  adapt_wavelet ((scalar*){f,u}, (double[]){f_tol,ux_tol,uy_tol},maxlevel , maxlevel-3 );
}

