@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
@include "_boundarydecl.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "couetteflow-channel.c"
// This code is written to implement a channel flow in the bottom  of the domain with width 2 and length 8. wiht the help of embed.h header and by giving the boundary conditions on u.n[embed.h](contraary to assuming u.n[] is the normal component. ).
// The code now has to be followed up to introduce an interface in the middle of the domain.  
//Author : Anvesh
//
//status: working. 
//comment: In this code giving the initial velocity as 0.0 makes it converge too quickly hence as a workaround the initial velocity is given as 0.001 and the mesh size has been refined to 10. 
//
// Date : 12-aug-2021 
//
//libraries used - 


#include "embed.h"
#include"navier-stokes/centered.h"
#include"vtk.h"


char name[100];
char name_vtk[100];
double H0, U0; 

double Reynolds = 20.;
int maxlevel = 7;
face vector muv[];

int main() 
{	       // Main program begins here
        L0 = 8.;            // Size of the square box



        H0 = 1.;            // Height of the channel
        U0 =10.;             // Velocity of the plate
        origin (-L0/2., -L0/8.);  // Origin is at the bottom centre of the box
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
u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left]    = dirichlet(0.);  // We give no pressure gradient to check for linear profile 
pf[left]   = dirichlet(0.);


u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow 
pf[right]  = dirichlet(0.);








event init (t = 0)
{
	boundary(all);

     foreach()
          
          u.x[] = 0.001 ;
  
     vertex scalar phi[];
  
     double eps = L0/(1 << 7)/1000.;
  
     foreach_vertex()
	     phi[] = intersection (-(y - L0/8. + eps), -(- L0/8. + eps - y));
	fractions (phi, cs, fs);
	boundary (all );

	mu = fm;




}


u.t[embed] = y > 0.5 ? dirichlet(0.) : dirichlet(0.);
u.n[embed] = y > 0.0 ? dirichlet(U0) : dirichlet(-U0);

// The above embed boundary condition is given as an IF condition to specify two different 
// embed conditions namely, the top embed wall going to the right and the
//  bottom embed wall going to the left. 

// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);






// Prost processing the results
event parview (i += 5  ; t <=5)
{
        scalar  m[];
        foreach()
                m[] = cs[]  ;
        boundary ({m});


        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);


}
// Using adaptive grid based on velocity
event adapt (i++) {
        adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 10);  //changes 
}


#endif
