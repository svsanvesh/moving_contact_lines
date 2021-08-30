

//This code simulates a couette flow in two concentric cylinders. 
//Written from : http://basilisk.fr/src/test/couette.c 
//Author: Anvesh 
// Status : Doesn't compile
//We need to link various libraries for *******view.h******** header.
//
//The method to compile the code is --   qcc -Wall -O2 [program].c  -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//Libraries included
//


#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"


int main()
{
  origin (-L0/2., -L0/2.);
  
  stokes = true;
  TOLERANCE = 1e-5;
  for (N = 16; N <= 256; N *= 2)
    run();
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {

	mu=fm ;
	 vertex scalar phi[];
	foreach_vertex()
	    phi[] = difference (sq(0.5) - sq(x) - sq(y), sq(0.25) - sq(x) - sq(y));
	boundary ({phi});
	fractions (phi, cs, fs);

	u.n[embed] = dirichlet (x*x + y*y > 0.14 ? 0. : - y);
	u.t[embed] = dirichlet (x*x + y*y > 0.14 ? 0. :   x);
	 foreach()
 	   un[] = u.y[];
}



event logfile (t += 0.01; i <= 1000) 
{
	  double du = change (u.y, un);
	  if (i > 0 && du < 1e-7)
	    return 1; /* stop */
}



#define powerlaw(r,N) (r*(pow(0.5/r, 2./N) - 1.)/(pow(0.5/0.25, 2./N) - 1.))

event profile (t = end)
{
  scalar utheta[], e[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (cs[] > 0.) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      e[] = utheta[] - powerlaw (r, 1.);
    }
    else
      e[] = p[] = utheta[] = nodata;
  }

  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  dump();
  
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("utheta", spread = -1);
  save ("utheta.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e.png");

  if (N == 32)
    foreach() {
      double theta = atan2(y, x), r = sqrt(x*x + y*y);
      fprintf (stdout, "%g %g %g %g %g %g %g\n",
	       r, theta, u.x[], u.y[], p[], utheta[], e[]);
    }
}
















