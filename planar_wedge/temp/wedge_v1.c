//This code is written to simulate a flow inside a wedge shaped 
//triangular region.
//Author : anvesh
// Method to compile code : qcc -O2 wedge_v1.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//Status :
// Domain is 1x1 box where we insert an embedded region.  
// Date : 11-oct-2021.
// libraries used- 


#include "embed.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "view.h" 


char name_vtk[100];
int maxlevel =6;
double Uplate;
double Reynolds = 2.0;       // Reynolds number
face vector muv[];		//viscosity. 
int main()
{

	L0=1;
	Uplate= 10.;
	origin (0., -L0);  // Origin is at thetop left corner. (middle ) 
	N=64; 
	mu= muv;
	run(); 



}
event properties (i++)
{
	foreach_face()
	 muv.x[] = Uplate*L0/Reynolds;
}

vertex scalar phi[];
event init (t = 0)
{

     foreach()

          u.x[] = 0.001 ;



     foreach_vertex()
             phi[] = -(y + x);
        boundary ({phi} );
        fractions (phi, cs, fs);

}
// Setting the boundary conditions
u.n[left] = dirichlet(0.0);
u.t[left] = dirichlet(-Uplate);


// outflow bc according to http://basilisk.fr/src/test/swirl.c
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);

// noslip bc on the inclined wall

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

//u.n[embed] = (y+x) < 0.0 ? neumann(0.) : dirichlet(0.);
//u.t[embed] = (y+x) < 0.0 ? neumann(0.) : dirichlet(0.);




// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);

// Post processing results  
event movies (i += 1  ; t <= 2)
{
	clear();
	  draw_vof ("phi");
	  cells();
//	  filled();
//	  isoline();
	  box();
	  save ("csscalar.mp4");
	  
	  save ("csscalar-.png");

	 clear();
	 vectors ("u");
	    box ();
	    save ("u.mp4");

scalar phi0[];
	foreach()
		phi0[]=phi[];

	foreach()
                sprintf (name_vtk, "data-%g.vtk", t);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,p,phi0},N,fpvtk,1);
}


// Using adaptive grid based on velocity
event adapt (i++) {
        adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 9); 
}

