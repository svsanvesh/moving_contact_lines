// THIS CODE IS WRITTEN TO REWRITE HE COUETTE WITH INTERFACE CODE AND RUN IT.
// AUTHOR: ANVESH 
// DATE : 9-SEPT-2021
// COMMETS: 
//FOR COMPILING : qcc new-couette-with-interface.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

//LIBRARIES USED- 

#include "navier-stokes/centered.h"
#include "vtk.h"
#include "vof.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"
#include "contact.h"
#include "tension.h"

// Computational parameters
double Reynolds = 20.0;       // Reynolds number
int maxlevel = 7;              // Maximum mesh refinement
face vector muv[];             // viscosity
double H0;
double U0;
char name_vtk[100];
double theta_bot = 90;
double theta_top = 90;
vector h[];
scalar f[], * interfaces = {f};
double rho1 = 1000., mu1 = 0.01, rho2 = 100., mu2 = 0.01;
h.t[bottom] = contact_angle (theta_bot*pi/180.);

int main() {                // Main program begins here
        L0 = .2;            // Size of the square box



//        H0 = 1.;            // Height of the channel
        U0 = 1.;             // Velocity of the bottom plate
        origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
        N = 128;
        mu = muv;           // constant viscosity. Exact value given below
	
	f.sigma = 1.;
        f.height = h;


        run();

}


/*
// Setting viscosity in the domain
event properties (i++)
{
        foreach_face()
         muv.x[] = U0*H0/Reynolds;
}
*/
// Setting the boundary conditions
u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left]    = dirichlet(0.);  // We give no pressure gradient to check for linear profile
pf[left]   = dirichlet(0.);


u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow
pf[right]  = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(-U0);


event init (t=0)
{
  fraction (f , x );
        boundary(all);



     foreach()
	     u.x[] =  0.01;

}




// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);

// Produce vorticity animation
event movies (i += 5  ; t <= 0.03)
{
	clear();
	  draw_vof ("f");
	  cells();
	  box();
	  save ("fscalar.mp4");


	 clear();
	 vectors ("u");
	    box ();
	    save ("u.mp4");






	foreach()
                sprintf (name_vtk, "data-%d.vtk", i);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,p,f},N,fpvtk,1);


}

/*

// Using adaptive grid based on interface position
event adapt (i++) {
        adapt_wavelet ((scalar*){f,u}, (double[]){0.1, 0.1,0.1}, 9);   
}
*/

event adapt(i++){
 adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){0.01,0.01, 0.01}, 10);
}
                    
