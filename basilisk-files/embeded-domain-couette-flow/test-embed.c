//This code is to test out the code in:  http://basilisk.fr/src/test/uf.c 
// code compilation : qcc test-embed.c -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//
//
//
//


#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

int main()
{
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  
  stokes = true;
  DT = 2e-5;
  TOLERANCE = HUGE;
  NITERMIN = 10;
  N = 32;

  run();
}

event init (t = 0)
{
  vertex scalar phi[];
  double eps = L0/(1 << 7)/1000.;
  foreach_vertex()
    phi[] = union (y - L0/4. + eps, - L0/4 + eps - y);
  boundary ({phi});
  fractions (phi, cs, fs);

  mu = fm;

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  foreach()
    u.y[] = 1.;
  boundary ((scalar *){u});
}

event logfile (i++; i <= 100)
{ 
	/*
  fprintf (stderr, "%d %d %d %d %d %d %d %.3g %.3g %.3g %.3g\n",
	   i,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resa*dt, mgu.resa, normf(u.y).max, normf(p).max);
*/
  foreach()
    if (x < -L0/2. + L0/N)
      printf (" %g\n", cs[]);
 // printf ("\n");
}

event profile (t = end)
{
  p.nodump = false;
  dump();
  assert (normf(u.y).max < 1e-3);
}
