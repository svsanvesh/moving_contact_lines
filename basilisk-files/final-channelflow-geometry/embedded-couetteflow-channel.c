// This code is written to implement a channel flow in the middle of the domain with width 2 and length 8. 
//Author : Anvesh 	
//status: Working (cp embeded-couetteflow-channel.c) 
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
int maxlevel = 6;
face vector muv[];

int main() 
{	       // Main program begins here
        L0 = 8.;            // Size of the square box



        H0 = 1.;            // Height of the channel
        U0 =10.;             // Velocity of the bottom plate
        origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
        N = 64;
        mu = muv;           // constant viscosity. Exact value given below

 U0 =10.;             // Velocity of the bottom plate
        origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box
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








event init (t = 0)
{
	boundary(all);

     foreach()
          
          u.x[] = cs[] ? 1.0  : 0.;
  
     vertex scalar phi[];
  
     double eps = L0/(1 << 7)/1000.;
  
     foreach_vertex()
	     phi[] = intersection (-(y - L0/8. + eps), -(- L0/8. + eps - y));
	fractions (phi, cs, fs);
	boundary (all );

	mu = fm;

	u.t[embed] = y > 0.5 ? dirichlet(0.) : dirichlet(0.);
	u.n[embed] = y > 0.0 ? dirichlet(U0) : dirichlet(-U0);
// The above embed boundary condition is given as an IF condition to specify two different 
// embed conditions namely, the top embed wall going to the right and the
//  bottom embed wall going to the left. 



}



// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);






// Produce vorticity animation
event parview (i += 10  ; t <=4)
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
        adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 7);  //channges 
}

