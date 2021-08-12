/* This file solves  for velocity and pressure field in a single phase Couette flow
 * The bottom plate is moving to the right and the to plate is moving to the left */
/* Date: 12-July-2021 */
// Author:  Harish(prime) Anvesh  
// Code comments: The code runs and provides the right results i.e.  the desired linear velocity profile.  
//
//THE VARIOUS CHANGES DONE IN  THE CODE AS AN EXPERIMENT-
//1. THE PRESSURE ON BOTH THE ENDS IS INTRODUCED
//2. THE LENGTH OF THE DOMAIN IS VARIED FOR VISUALIZATION IN THE FORM OF L0. 
//3. THE event MOVIES FUNCTION PRINITNG TIMESTEP IS INCREASED TO INCREASE THE EXECTION SPEED. 
//4. THE PLATE VELOCITY IS SIGNIFICANTLY INCREASED OT SEE THE VISUIEAL EFFECTS -  {10, 5, 2.5  }
//5. THE MASK FUNCTION m[] = cs  is modified to see the changes but wthout without any success 
//6. IN EVENT INIT(INITIALIZATION STEP) THE U.X[] IS HARDCODED TO BE 0 EVERYWHERE.  (this crashes the code and does not compile)           
//7. THE PRESSURE IS SET TO (5.) ON BOTH LEFT AND RIGHT SIDE. 
//8. THE TEST CASE IS RUN FOR A LONGER TIME IN THE HOPE TO GET STEADY STATE LINEAR VELOCITY PROFILE.(IN PREVIOUS CASE THE BACKWARD MOVEMENT  OF VELOCITY WAS STABILISING)[T=30] 
//9. THE GRID RESOLUTION IS SIGNIFICANTLY INCREASED SO AS TO AVAID THE SEGEMENTATION ERROR THAT WE SAW PREVIOUSLY [MAXLEVEL= 15, N=32768] {too much doesn't work} 
//10. THE GRID RESOLUTION IS SIGNIFICANTLY INCREASED SO AS TO AVAID THE SEGEMENTATION ERROR THAT WE SAW PREVIOUSLY [MAXLEVEL= 13, N=8192]  
//11. THE INITIAL VELOCITY EARLIER SPECIFIED IS REMOVED TO CHECK ITS EFFECT. [IT GIVES AS ERROR. ] 
//12. WE GIVE NEUMANN CONDITIONS FOR PRESSURE SO AS TO MAKE THE FLOW CONTINUOUS. 
//13. WE MOVE THE TOP PLATE WITH THE SAME VEOCITY TO MAKE A COUETTE FLOW WITH ZERO RATE OF CHANGE OF MASS. 
//14. IN THE BOUNDARY CONDITIONS THE TOP PLATE IS MOVED WOTH THE HELP OF THE [TOP] COMMAND ITSELF----THE [EMBED] DOES NOT WORK.  
//
//
// CONCLUSION: tHE CODE WORKS FOR A 1X1 DOMAIN WITH BOUNDARY CONDITIONS GIVEN AS SO AND WE ARE ABLE TO PLOT A LINEAR VELOCITY PROFILE. 
//
//
//
//
//// Libraries included
#include "embed.h"
#include "navier-stokes/centered.h"
#include "vtk.h"

// Computational parameters
double Reynolds = 100.0;       // Reynolds number
int maxlevel = 5;              // Maximum mesh refinement
face vector muv[];             // viscosity
double H0;
double U0;
char name[100];
char name_vtk[100];

//#define Rectangle ()

int main() {                // Main program begins here
	L0 = 1.;            // Size of the square box
	


	H0 = 1.;            // Height of the channel
	U0 =10;             // Velocity of the bottom plate
	origin (-L0/2, 0.0);  // Origin is at the bottom centre of the box
	N = 32 ; 
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
//	vertex scalar phi[];
//	foreach_vertex() {
//phi[] = intersection (0.5 - y, 0.5 + y);
//	 phi[] = (y < H0); 
//	}
//	boundary ({phi});
//	fractions (phi, cs, fs);

	boundary(all);

     foreach()

	  u.x[] = cs[] ? H0  : 0.;	//(IT DOESNT ) 	THIS LINE IS COMMENTED TO CHECK IF THE CODE STILL RUNS 

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

        sprintf (name, "vort-%g.ppm", t);
        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        FILE * fp = fopen (name, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);
        output_ppm (u.x, fp, min = -2, max = 2, n = 512);


}

// Using adaptive grid based on velocity
event adapt (i++) {
	adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 5);
 
