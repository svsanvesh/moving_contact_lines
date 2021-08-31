//THE CODE IMPLEMENTS ATOMISATION.C FROM : http://www.basilisk.fr/src/examples/atomisation.c 
//WRITTEN BY : ANVESH 
//THE CODE IS PRIMARILY WRITTEN TO CHECK THE ADAPT_WAVELET FUNCTION WHEN ONLY F0 IS GIVEN AND SOME FUNCTIONALITIES LIKE THE DROP COUNTING HAS BEEN COMMENTED. 
//DATE: 31-AUG-2021
//COMMENTS : IT IS WORKING 
//STATUS : WORKING . 
//
//
//
//LIBRARIES USED - 



#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

#define radius 1./12.
#define length 0.025
#define Re 5800
#define SIGMA 3e-5

int maxlevel = 9;
double uemax = 0.1;

scalar f0[];
u.n[left]  = dirichlet(f0[]*(1. + 0.05*sin (10.*2.*pi*t)));
u.t[left]  = dirichlet(0);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = f0[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  init_grid (64);
  origin (0, -1.5, -1.5);
  size (3.);

  rho1 = 1., rho2 = 1./27.84;
  mu1 = 2.*radius/Re*rho1, mu2 = 2.*radius/Re*rho2;  
  f.sigma = SIGMA;

  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {

    refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel);

    fraction (f0, sq(radius) - sq(y) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels

    foreach() {
      f[] = f0[]*(x < length);
      u.x[] = f[];
    }
    boundary ({f,u.x});
  }
}

event logfile (i++) {
  if (i == 0)
    fprintf (stderr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (stderr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
}

event movie (t += 1e-2)
{
#if dimension == 2
  scalar omega[];
  vorticity (u, omega);
  view (tx = -0.5);
  clear();
  draw_vof ("f");
  squares ("omega", linear = true, spread = 10);
  box ();
#else // 3D
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid}); // not used for the moment
  view (camera = "iso",
	fov = 14.5, tx = -0.418, ty = 0.288,
	width = 1600, height = 1200);
  clear();
  draw_vof ("f");
#endif // 3D
  save ("movie.mp4");
}

event snapshot (t = 0.1; t += 0.1; t <= 3.8) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid});
  dump (name);
}


event adapt (i++) {
  adapt_wavelet ({f}, (double[]){0.01}, maxlevel);
}

