//This is a simualtion to visualize the flow field near a moving contact line. 
//The geometry of the problem is a sqaure domain of size L=5*l_c ; where l_c = 3. (ALL LENGTHS IN mm)
//it is a 15x15 square with interface in the middle, horizontally. 
//Author- Anvesh 
//Date - 15-Nov-2021
//Comments: 
//Status : working 
//
//
//Libraries used - 

//#include "navier-stokes/conserving.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "adapt_wavelet_leave_interface.h"
#include "contact.h"
#include "tension.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"

double Reynolds = 2.0;       // Reynolds number
int maxlevel = 9;              // Maximum mesh refinement
//face vector muv[];             // viscosity
char name_vtk[100];		// vtk file name decleration.
double U0;
double H0;


/*
	#define grav  9.81 // gravitational acceleration
        #define rhoL 1000  //density of water
        #define muL 0.001 //viscosity of water
        #define surf 0.072 // surface tension air-water
        #define rhoG 1.225 //density of air
        #define muG 0.0000181 // viscosity of air

*/

vector h[];  //HEIGHT FUNCTION 
double theta0 ; 
int main() 
{	
        L0 = 15.;            // Size of the square box

//        H0 = 1.;            // Height of the channel
	dt=0.1;
        U0 = 1 ;             // Velocity of the left plate
        origin (0, -L0/2);  // Origin is at the bottom centre of the box
        N = 256;
      //  mu = muv;           // constant viscosity. Exact value given below

	stokes = true;
        f.sigma = 1.;
        f.height = h;
	display_control (maxlevel, 6, 12);

	theta0 = 150; 
	h.t[left] = contact_angle (theta0*pi/180.); // Left contact angle near the moving wall 
	h.t[right] = contact_angle (pi/2);  // right contact angle of 90 degrees. 
	
	//The viscosity and desinties of the two fluids is specified here. 
	rho1 = rhoG;
	mu1 = muG;
	rho2 = rhoL;
	mu2 = muL;
	
	/*
	const face vector muc[] = {.1,.1};
	mu = muc;
	*/
        run();

}

event adapt (i++) {
scalar impose_refine[], f1[];
foreach(){
          //        f1[] = f[];
                  if(y>-1.&&y<1.&&x<15)
                          impose_refine[] = noise();
                  else
                          impose_refine[] = 0;
}  
          boundary({impose_refine});
 adapt_wavelet2((scalar *){impose_refine},(double []){0.001},(int []){13,10},3);

}
event init (t = 0)
{
	fraction (f , y );
	boundary ({f});
     
  foreach()
             u.x[] =  0.0001;
}



// gravity is given in the vertically down direction.(-9.81)
event acceleration (i++)
{
        face vector av = a;
        foreach_face(x)
                av.y[] = -9.81;
}


// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(-U0);


u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.0);



// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);




// Produce vorticity animation
event movies (i += 5  ; t <= 20.)
{
        dump( );
        foreach()
                sprintf (name_vtk, "data-%d.vtk", i);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,p,f},N,fpvtk,1);

}
event adapt (i++) {
        adapt_wavelet ((scalar*){f,u}, (double[]){0.1, 0.1,0.1}, 9,6);
}
/*
// Using adaptive grid based on interface position
event adapt(i++){
 adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){0.0 ,0.0, 0.01}, maxlevel);
}

*/
