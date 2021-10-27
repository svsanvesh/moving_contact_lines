// THIS IS A LID DRIVEN CAVITY CODE FOR CHECKING THE DEFAULT PRESSURE BOUNDARY CONDITION TAKEN BY 
// BASILISK. 
// AUTHOR : ANVESH 
// STATUS : WOKRING 
//
// GEOMETRY -  FOUR WALLS LEFT WALL MOVING  DOWN. 
//






#include "navier-stokes/centered.h"
#include "vtk.h"
#include "view.h"

#define rho  1000; 

char name_vtk[100];
int maxlevel =6;
double U0 ; 
double Reynolds = 2.0;       // Reynolds number
face vector muv[];              //viscosity.
int main()
{
	U0=1.0;
        L0=1. ;
        origin (0., -L0);  // Origin is at thetop left corner.
        N=256;
        mu= muv;
        run();

}


event properties (i++)
{
        foreach_face()
         muv.x[] =U0*L0/Reynolds;
}

event init (t = 0)
{

     foreach()

          u.x[] = 0.001 ;


}

// left
	u.t[left]   = dirichlet(-U0);

	p[left]	   =  neumann(0.);	
// top
	u.t[top]    = dirichlet(0.); 
	p[top]	   =  dirichlet(1.);	
//right
	u.t[right]  = dirichlet(0.);

	p[right]   =  neumann(0.);	

//bottom 
	u.t[bottom] = dirichlet(0.);

	p[bottom]   =  neumann(0.);	


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


