/* This file solves  for velocity and pressure field in a single phase Couette flow
 * The bottom plate is moving to the right and the to plate is moving to the left */
/* Date: 12-July-2021 */
// Author: Anvesh  
// Code comments: The code runs and provides the right results i.e.  the desired linear velocity profile.  
// CONCLUSION: THE CODE WORKS FOR A 1X1 DOMAIN WITH BOUNDARY CONDITIONS GIVEN AS SO AND WE ARE ABLE TO PLOT A LINEAR VELOCITY PROFILE. 
//
//
//
//
//// Libraries included
//#include "embed.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include"two-phase.h"


// Computational parameters
double Reynolds = 20.0;       // Reynolds number
int maxlevel = 5;              // Maximum mesh refinement
face vector muv[];             // viscosity
double H0;
double U0;
char name_vtk[100];
scalar f[], * interfaces = {f};
//#define Rectangle ()

int main() {                // Main program begins here
	L0 = 1.;            // Size of the square box
	


	H0 = 1.;            // Height of the channel
	U0 =1;             // Velocity of the bottom plate
	origin (-L0/2., -L0/2.0);  // Origin is at the bottom centre of the box
	N = 32; 
	mu = muv;           // constant viscosity. Exact value given below

	run();

}

// Setting viscosity in the domain
event properties (i++)
{
	foreach_face()
	 muv.x[] = U0*H0/Reynolds;
	//muv.x[] = fm.x[]/Reynolds; 
}

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
u.t[top] = dirichlet(U0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(-U0);


event init (t=0)
{
	fraction (f , x );
	boundary(all);

     foreach()

	  u.x[] = cs[] ? 1.0 : 0.;	//(IT DOESNT ) 	THIS LINE IS COMMENTED TO CHECK IF THE CODE STILL RUNS 

}

// Printing out standard text outputs on the screen
event logfile (i++)
	fprintf (stderr, "%d %g\n", i, t);

// Produce vorticity animation
event movies (i += 10  ; t <=10)
{
	scalar omega[], m[];
	vorticity (u, omega);
	foreach()
		m[] = cs[]  ;
	boundary ({m});
			

        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        output_vtk ({u.x,u.y,p,f},N,fpvtk,1);


}

// Using adaptive grid based on velocity
event adapt (i++) {
	adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 10);  //channges 
}
