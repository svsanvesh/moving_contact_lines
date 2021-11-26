//This is a simualtion to visualize the flow field near a moving contact line. 
//The geometry of the problem is a sqaure domain of size L=5*l_c ; where l_c = 3. (ALL LENGTHS IN mm)
//it is a 15x15 square with interface in the middle, horizontally. 
//Author- Anvesh 
//Date - 13-Nov-2021
//Comments: There is a problem with draining of the bottom fluid.  
//Status :  not working 
//
//
//Libraries used - 

//#include "navier-stokes/conserving.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "adapt_wavelet_leave_interface.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"



double Reynolds = 2.0;       // Reynolds number
int maxlevel = 9;              // Maximum mesh refinement
face vector muv[];             // viscosity
char name_vtk[100];		// vtk file name decleration.
scalar f[], * interfaces = {f};
double rho1 = 1., mu1 =1.0,rho2 = 10., mu2 = 50.0;
double U0;
double H0;
        #define grav  9.81 // gravitational acceleration
	#define rhoL 1000  //density of water
        #define muL 0.001 //viscosity of water
        #define surf 0.072 // surface tension air-water
        #define rhoG 1.225 //density of air
        #define muG 0.0000181 // viscosity of air


vector h[];  //HEIGHT FUNCTION 
double theta0 = 150; 
h.t[left] = contact_angle (theta0*pi/180.); // Left contact angle near the moving wall 
h.t[right] = contact_angle (pi/2);  // right contact angle of 90 degrees. 

int main() {                // Main program begins here
        L0 = 15.;            // Size of the square box

	CFL = 0.45;


//        H0 = 1.;            // Height of the channel
        U0 = 1.;             // Velocity of the bottom plate
        origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
        N = 256;
//        mu = muv;           // constant viscosity. Exact value given below



  	stokes = true;
        f.sigma = surf;  // surface tension for the interface is given here. 
        f.height = h;
	display_control (maxlevel, 6, 12);


	const face vector muc[] = {.1,.1};
	mu = muc;
        run();

}

// gravity is given in the vertically down direction.(-9.81)
event acceleration (i++)
{
	face vector av = a;
	foreach_face(x)
		av.y[] = -grav;
}
event init (t = 0)

{
	fraction (f , y );
	boundary ({f});
     
  foreach()
             u.x[] =  0.0001;
}


// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(-U0);
p[left]    = dirichlet(0.);  // We give no pressure gradient to check for linear profile
pf[left]   = dirichlet(0.);


u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow
pf[right]  = dirichlet(0.);



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
        adapt_wavelet ((scalar*){f,u}, (double[]){0.1, 0.1,0.1}, 10,7);
}

// Using adaptive grid based on interface position
event adapt(i++){
 adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){0.01,0.01, 0.01}, maxlevel);
}

