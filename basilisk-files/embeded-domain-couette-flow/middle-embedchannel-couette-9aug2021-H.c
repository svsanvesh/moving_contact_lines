// This code is written to visualize the fluid flow domain and to find a way to change the fluid domain from square to rectangular.And this is done with the help of phi[] fractions and then taking intersection.  
// author : Anvesh & Harish 
// Status : working
// This is based on  Von-karman.c code given in http://basilisk.fr/src/examples/karman.c as a base on top of which changes are made.  






#include "embed.h"
#include"fractions.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"
#include"vtk.h"



char name[100];
char name_vtk[100];


int  H0= 2. ;



scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 5;
face vector muv[];
int main() {
  L0 = 8.;
  origin (-L0/2, -L0/2);
  N = 32;
  mu = muv;
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);

  run();
}
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/Reynolds;
}
u.n[left]  = neumann(0.);
u.t[left]  = neumann(0.);
p[left]    = dirichlet(0.);
pf[left]   = dirichlet(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

//Dirichlet couditions-
u.n[embed] = y > 0.25 ? dirichlet(0.):  dirichlet(0.);
u.n[embed] = y < -0.25 ? dirichlet(0.):  dirichlet(0.);
u.t[embed] = y > 0.25 ? dirichlet(2.):  dirichlet(0.);
u.t[embed] = y < -0.25 ? dirichlet(-2.):  dirichlet(0.);



  vertex scalar phi[];
// Initial condition
event init (t = 0)
{

  //vertex scalar phi[];
      foreach_vertex(){
	      if(fabs(y) <= 1.0)
	                phi[] = 1.0;
		
              else
                        phi[] = 0.0;
        }

      foreach(){
	      if(fabs(y) <= 1.0)
		         u.x[] = y;
              else
		        u.x[] = 0.0;
        }


//  foreach_vertex() {
//    phi[] = intersection (0.5 - y, 0.5 + y);
//  }
  boundary ({phi});
  fractions (phi, cs, fs);

//  foreach()
 //   u.x[] = 0.;
}

event logfile (i++)
{
       //	fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
 	foreach()
	       printf("%g \n  ",cs[])	
} 
  
event movies(i += 4; t <= 20.)
{
  scalar omega[], m[], phi0[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
         foreach()
          phi0[] = phi[];


  //	sprintf (name, "vort-%g.ppm", t);
        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
    //    FILE * fp = fopen (name, "w");
        output_vtk ({u.x,u.y,p,phi0},N,fpvtk,1);
      //  output_ppm (u.x, fp,min= -10 , max =10, n = 512);



}
event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}
