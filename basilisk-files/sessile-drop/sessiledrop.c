//THIS CODE IS WRITTEN TO REPLICATE THE CODE OF SESSILE.C FROM :http://basilisk.fr/src/test/sessile.c 
//





#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};

vector h[];
double theta0 = 30;
h.t[bottom] = contact_angle (theta0*pi/180.);

int main()
{
  size (2);

  const face vector muc[] = {.1,.1};
  mu = muc;

  f.height = h;

  f.sigma = 1.;

  for (theta0 = 15; theta0 <= 165; theta0 += 15)
    run();
}

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

event end (t = 10)
{
  output_facets (f, stdout);
  
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}
