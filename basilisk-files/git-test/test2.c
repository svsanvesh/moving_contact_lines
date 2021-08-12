//
//
//
//qcc -O2 test.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm




#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

vector h[];

#define MAX_LEVEL 7
#define TIME_STEP 1./25 // Timestep for movie and log file.
#define MAX_TIME 10.0 // Time to stop at.
#define beta 0.93 // Location of the fluid interface

// Fluid 1 is water at 25°C
#define rhoWater 0.9970
#define muWater 0.8937
// Source: retrieved 7/9/20 from
// https://en.wikipedia.org/wiki/Water_(data_page)#Liquid_physical_properties

// Fluid 2 is sunflower oil at 25°C
#define rhoOil 0.9188
#define muOil 49.14
// Source: retrieved 7/9/20 from
// https://en.wikipedia.org/wiki/Sunflower_oil#Physical_properties

#define Fr 1./sqrt(9.81)
#define We 1.*1.*1./0.037

int main () {
  L0 = 8.;
  origin(-L0/2., -0.5);
  init_grid(1 << MAX_LEVEL);
  // Set up fluids, doing this in the init event doesn't work.

  // Fluid 1 is water 
  rho1 = rhoWater;
  mu1 = muWater;

  // Fluid 2 is sunflower oil
  rho2 = rhoOil;
  mu2 = muOil;

  // Surface tension coefficient
  f.sigma = 1./We;
  f.height = h;
  run();
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1./sq(Fr);
}

event init(t = 0) {

  // Flow out of the right -> into the left.
  periodic(right);

  // positive -> fluid 1
  // negative -> fluid 2 
  // bottom beta% is oil.
  fraction (f, beta-y + 0.04*sin(4*pi*x));

  // Initially velocity is 0 everywhere.
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }

  mask (y > 0.5 ? top : none);
  mask (y < -0.5 ? bottom : none);
  
  // No-slip boundary conditions.
  u.t[top] = dirichlet(1.);
  u.t[bottom] = dirichlet(0.);

}

event vof (i++, first);

event amplitude (t += 1) {
  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t, max);
  fflush (fp);
  //se += sq(max - prosperetti[ne][1]); ne++;
}



event adapt (i++) {
  // error thresholds for            f     u.x   u.y
  adapt_wavelet ({f, u}, (double []){1e-2, 1e-2, 1e-2}, MAX_LEVEL);
}

event profiles (t = 0; t+=1.0; t<=1000) // RC restricted the output a little, don't overdo it at first!
{
  FILE * fp = fopen("xprof", "a");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g %g\n", t, y, interpolate (u.x, 0, y));
  fclose (fp);

  fp = fopen("yprof", "a");
  for (double x = -4; x <= 4; x += 0.01)
    fprintf (fp, "%g %g %g\n", t, x, interpolate (u.y, x, 0));
  fclose (fp);

  fp = fopen('heights', "a");
  fprintf(fp, "%g %g\n", t, statsf(h).max);
  fclose (fp);
}
