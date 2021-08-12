/***** THIS CODE IS USEED TO IMPLEMENT EMBED BOUNDARY CONDITIONS OF DIRICHLET MOVING TOP PLATE ******/

/* This file solves
 * for velocity and pressure field in a single phase Couette flow
 * The bottom or embed(top) plate moving */
/* Date: 5-August-2021 */
// Auther : Anvesh 
// status: working
// Code comments : Compiling and working code(with success in embed condtions )
//// Libraries included

#include "embed.h"
#include "navier-stokes/centered.h"
#include "vtk.h"

// Computational parameters
double Reynolds = 20.0;
int H0; 
int maxlevel = 5;              // Maximum mesh refinement
face vector muv[];             // viscosity
double U0;
char name[100];
char name_vtk[100];

//#define Rectangle ()

int main() {                // Main program begins here
	L0 = 8.;            // Size of the square box
	H0 = 1.;            // Height of the channel
	U0 = 10.  ;            // Velocity of the bottom plate
	origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
	N = 32 ; 
	mu = muv;           // constant viscosity. Exact value given below

	run();
}

// Setting viscosity in the domain
event properties (i++)
{
	foreach_face()
	 muv.x[] = U0*H0/Reynolds;
}


vertex scalar phi[];

event init (t=0)
{

  foreach_vertex() {
    phi[] = intersection (y + 0.5, y - 0.5);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);


     foreach()
	  u.x[] = 0.;	// Zero velocity field inside the channel


}


// Setting the boundary conditions
u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left]   = dirichlet(0.);  // We give no pressure gradient to check for linear profile 
pf[left]  = dirichlet(0.0);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow 
pf[right]  = dirichlet(0.);


//Dirichlet couditions- 
u.n[embed] = y > 0.5 ? dirichlet(0.):  dirichlet(0.);
u.n[embed] = y < -0.5 ? dirichlet(0.):  dirichlet(0.);
u.t[embed] = y > 0.5 ? dirichlet(U0):  dirichlet(0.);
u.t[embed] = y < -0.5 ? dirichlet(-U0):  dirichlet(0.);


// Printing out standard text outputs on the screen
event logfile (i++)
	fprintf (stderr, "%d %g \n", i, t);

// Produce vorticity animation
event movies (i += 5  ; t <=3)
{
	scalar omega[], m[], phi0[];
	vorticity (u, omega);
	foreach()
		m[] = cs[];
	boundary ({m});

        foreach()
          phi0[] = phi[];


//        sprintf (name, "vort-%g.ppm", t);
        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        FILE * fp = fopen (name, "w");
        output_vtk ({u.x,u.y,p,phi0},N,fpvtk,1);

	
}

// Using adaptive grid based on velocity
event adapt (i++) {
	adapt_wavelet ({cs,u}, (double[]){1e-3,3e-3,3e-3}, maxlevel, 6);
}


