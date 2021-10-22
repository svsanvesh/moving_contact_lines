// THIS CODE IS WRITTEN TO SIMULATE A FLOW INA WEDGE WITH LEFT FACE MOVING DOWN WITH THE CENTRE BEING IN THE TOP LEFT CORNER. 
//
//
//Author : Anvesh
//
//status: working. 
//comment: In this code giving the initial velocity as 0.0 makes it converge too quickly hence as a workaround the initial velocity is given as 0.001 and the mesh size has been refined to 10. 
//
// Date : 21-oct-2021 
//
//libraries used - 


#include "embed.h"
#include"navier-stokes/centered.h"
#include"vtk.h"


char name[100];
char name_vtk[100];
double H0, U0; 

double Reynolds = 0.2;
int maxlevel = 7;
face vector muv[];

int main() 
{	       // Main program begins here
        L0 = 1.;            // Size of the square box



        H0 = 1.;            // Height of the channel
        U0 =1.;             // Velocity of the bottom plate
        origin (0, -L0/2);  // Origin is at the bottom centre of the box
        N = 128;
        mu = muv;           // constant viscosity. Exact value given below

        run();


}
// Setting viscosity in the domain

event properties (i++)
{
        foreach_face()
         muv.x[] = U0*H0/Reynolds;
}

// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);


// Outflow BC according to http://basilisk.fr/src/test/swirl.c
u.n[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);








event init (t = 0)
{

     foreach()
	     u.x[] = 0.001 ;
  
     vertex scalar phi[];
  
  
     foreach_vertex()
	     phi[] = -(y+x) ;
	fractions (phi, cs, fs);
	boundary (all );

}


u.t[embed] =  dirichlet(0.);
u.n[embed] =  dirichlet(0.0);


// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);






// Prost processing the results
event parview (i += 10  ; t <=5)
{
        scalar  m[];
        foreach()
                m[] = cs[]  ;
        boundary ({m});


        sprintf (name_vtk, "data-%d.vtk",i);
        FILE * fpvtk = fopen (name_vtk, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);

	dump(); 
}
// Using adaptive grid based on velocity
event adapt (i++) {
        adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 10);  //changes 
}

