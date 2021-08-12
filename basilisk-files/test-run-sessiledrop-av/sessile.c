#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"
#include "vtk.h"
#include "embed.h"

scalar f[], * interfaces = {f} ;




char name[100];
char name_vtk[100];



vector h[];
double theta0 = 30;
h.t[bottom] = contact_angle (theta0*pi/180.);

int main()
{
  size (2);

//We use a constant viscosity.
  const face vector muc[] = {.1,.1};
  mu = muc;

// the height function is associacted with the VOF tracer. 
 f.height = h;

//We set the surface tension coefficient and run for the range of contact angles.
  f.sigma = 1.;

//  for (theta0 = 30; theta0 <= 150; theta0 += 30)
    run();
}

//The initial drop is a quarter of a circle.
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


event movies (i += 10 )
{
        scalar omega[], m[];
        vorticity (u, omega);
        foreach()
                m[] = cs[]  ;
        boundary ({m});
output_ppm (omega, file = "vort.mp4", box = {{-L0/2., 0.0},{L0/2., 1. }},

                min = -0.5, max = 2.0, linear = true, mask = m);

        sprintf (name, "vort-%g.ppm", t);
        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        FILE * fp = fopen (name, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);
        output_ppm (u.x, fp, min = -2, max = 2, n = 512);


}














event end (t = 10)
{
  output_facets (f, stdout);
  
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}


