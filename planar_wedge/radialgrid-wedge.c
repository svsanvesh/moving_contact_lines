//This code is written to simulate a flow inside a wedge shaped 
//triangular region.
//Author : anvesh
//THIS IS THE CASE WE TAKE ALL WALLS AS CLOSED AND NO SLIP IS APPLIED ON EACH WALL. 
// Method to compile code : qcc -O2 wedge-driven-cavity.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//Status :working 
// Domain is 1x1 box where we insert an embedded region.  
// Date : 20-oct-2021.
// libraries used- 
// Version - 5
// problem : 2


#include "embed.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "view.h" 
#include "radial.h"




int maxlevel =9;
double Uplate=1.0;
double Reynolds = 0.2;       // Reynolds number
face vector muv[];		//viscosity. 
int main()
{
     // maximum timestep
        DT = 0.1;
	L0=1;
	origin (-L0/2, -L0/2);  // Origin is at thetop left corner.
	N=256; 
	mu= muv;
	run(); 


}
event properties (i++)
{
	foreach_face()
        muv.x[] = Uplate*L0/Reynolds;
}
event adapt (i++) {
	scalar s[] ;   // Corner refinement scheme taken from : https://groups.google.com/g/basilisk-fr/c/EMg6USbSVq0/m/R79XWddDBAAJ 
	foreach()
	//	s[]= -(y+x); 
		s[] = sqrt(x*x+y*y)-sqrt(0.5);
		boundary ({s});
		adapt_wavelet ((scalar*){s}, (double[]){0.1},  11, maxlevel);
}
// Setting the boundary conditions
// Boundary conditions on moving plate
u.n[left] = dirichlet(0.0);
u.t[left] = dirichlet(-1.);

/* // Outflow BC according to http://basilisk.fr/src/test/swirl.c
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);
*/
// Bottom wall boundary conditions - slip wall
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);


 scalar phi[];   // We solve for flow inside a wedge
// Setting the initial conditions
event init (t = 0)
{
double eps = L0/(1 << 7)/10.;

	foreach_face() {
	     phi[] =-(y + x - eps );
     }
     boundary ({phi} );
     fractions (phi, cs, fs);

     foreach() {
   	 u.x[] = cs[] ? 0.001 : 0.;
     }

}

// No-slip BC on the inclined wall
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
p[embed]= dirichlet(1.);

event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);


char name_vtk[100];
char name[80];
// Post processing results  
event movies (i += 50  ; t <= 10.)
{ 


	scalar omega[], m[];
	vorticity (u, omega);
	boundary ({cs});
	output_ppm (omega, file = "vort.mp4", min = -2, max = 2, linear = true, mask = cs);
	
	
	sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);
	dump();

}

// Using adaptive grid based on velocity
event adapt (i++) {
        adapt_wavelet ((scalar*){u,phi}, (double[]){3e-2,3e-2,0.001}, 9,7); 
}
