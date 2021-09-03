/* This file solves  for velocity and pressure field in a single phase Couette flow
 * The bottom plate is moving to the right and the to plate is moving to the left */
/* Date: 31-Aug-2021 */
// Author: Anvesh  
// Code comments: The code runs and provides the right results i.e.  the desired linear velocity profile.  
// CONCLUSION: THE CODE WORKS FOR A 1X1 DOMAIN WITH BOUNDARY CONDITIONS GIVEN AS SO AND WE ARE ABLE TO PLOT A LINEAR VELOCITY PROFILE. 
//
// for compiling : qcc -O2 couette-with-interface.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//
//// Libraries included

#include "navier-stokes/centered.h"
#include "vtk.h"
#include "two-phase.h"
//#include "contact.h"
//#include "tension.h"
#include "view.h"

// Computational parameters
double Reynolds = 5.0;       // Reynolds number
int maxlevel = 6;              // Maximum mesh refinement
face vector muv[];             // viscosity
double H0;
double U0;
char name_vtk[100];

/*
double theta_bot = 90;
double theta_top = 90;
vector h[];
h.t[bottom] = contact_angle (theta_bot*pi/180.);
h.t[top] = contact_angle (theta_top*pi/180.);
*/


int main() {                // Main program begins here
	L0 = 1.;            // Size of the square box
	

	H0 = 0.2;            // Height of the channel
	U0 =10.0;             // Velocity of the bottom plate
	origin (-L0/2., -L0/2.0);  // Origin is at the bottom centre of the box
	N = 64;


/*        f.sigma = 1.;
	f.height = h;
*/
	run();

}

// Setting viscosity in the domain
event properties (i++)
{
	foreach_face()
	 muv.x[] = U0*H0/Reynolds;
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
u.t[top] = dirichlet(0.0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(-U0);


event init (t=0)
{
	fraction (f , x );
	boundary(all);

     foreach()

	  u.x[] = 0.01;	//(IT DOESNT ) 	THIS LINE IS COMMENTED TO CHECK IF THE CODE STILL RUNS 

}

// Printing out standard text outputs on the screen
event logfile (i++)
	fprintf (stderr, "%d %g\n", i, t);


//THIS IS WHERE GFS FILES AARE GENERATED. 
/*
event snapshot (i += 1  ; t <=0.4) {
  char name[80];
  sprintf (name, "snapshot-%d.gfs", i);
  scalar pid[];
  foreach()
    pid[] = tid();
  output_gfs (file = name);
  dump("t"); 
}
*/
event mp(i += 1  ; t <=0.5) 
{

  clear();
  draw_vof ("f");
  begin_mirror ({0,-1});
  cells();
  end_mirror();

  save ("movie.mp4");
}	

	
// Produce paraviewable files 
event movies (i += 1  ; t <=0.4)
{
	foreach()
	       	sprintf (name_vtk, "data-%d.vtk", i);
        	FILE * fpvtk = fopen (name_vtk, "w");
        	output_vtk ({u.x,u.y,p,f},N,fpvtk,1);


}
/*event snapshot (i += 1  ; t <=0.2) {
  char name[80];
  dump (name, list = {f});
} */


// Using adaptive grid based on interface position
event adapt (i++) {
	adapt_wavelet ((scalar*){f,u}, (double[]){0.1, 0.1,0.1}, 8);   
}







