// This code is written to simulate a wedge flow with right  wall moving down. 
// author : aNVESH 
//
// LIBRARIES USED - 
//
#include "navier-stokes/centered.h"
#include "embed.h"
#include "vtk.h"

face vector muv[];              //viscosity. 
double Reynolds = 0.2;       // Reynolds number

int U0=1;
// viscosity is defined here.
event properties (i++)
{
        foreach_face()
        muv.x[] = U0*L0/Reynolds;
}


// The main function is here 
int main()
{
        L0=1;
//        origin (-L0/2, -L0/2);  // Origin is at thetop left corner.
        N=128;
        mu= muv;
        run();
}

// Initial condition for the study. 
event init (t = 0 )
{
vertex scalar phi[];   // We solve for flow inside a wedge
	foreach()
		phi[]= -(y-x);


	boundary ({phi} );
	fractions (phi, cs, fs);

	foreach() 
	{
         u.x[] = cs[] ? 0.001 : 0.;
        }

}
// Boundary condition 
// right boundary 

u.n[right]= dirichlet(0.0);
u.t[right]= dirichlet(0.0);

// bottom boundary 
// Outflow BC according to http://basilisk.fr/src/test/swirl.c
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);

// embed boundary condition on the line y=x . 

u.n[embed]= dirichlet(0.0);
u.t[embed]= dirichlet(-U0);

char name_vtk[100];
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);

// Post processing results  
event movies (i += 1  ; t <= 0.8)
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









