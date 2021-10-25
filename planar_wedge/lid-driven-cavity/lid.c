#include "navier-stokes/centered.h"
#include "vtk.h"
#include "view.h"


char name_vtk[100];
int maxlevel =6;
double U0;
double Reynolds = 2.0;       // Reynolds number
face vector muv[];              //viscosity.
int main()
{
	U0=1.0;
        L0=1. ;
        origin (0., -L0);  // Origin is at thetop left corner. (middle )
        N=256;
        mu= muv;
        run();

}


event properties (i++)
{
        foreach_face()
         muv.x[] = U0*L0/Reynolds;
}

event init (t = 0)
{

     foreach()

          u.x[] = 0.001 ;


}
	u.t[left]   = dirichlet(-1.);
	u.t[top] = neumann(0.); 
	p[top] = dirichlet(1.);
	u.t[right]  = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);



// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);


// Post processing results  
event movies (i += 1  ; t <= 10)
{

scalar phi0[];

        foreach()
                sprintf (name_vtk, "data-%g.vtk", t);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,p},N,fpvtk,1);
		dump();
}
// Using adaptive grid based on velocity
event adapt (i++) {
        adapt_wavelet ((scalar*) {u}, (double[]){3e-2,3e-2}, 9, 7);
}


