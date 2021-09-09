#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"
#include "navier-stokes/perfs.h"

int minlevel = 4;
int maxlevel = 10;

double f_ = 0.75,
       s_ = 0.5334;
double uemax = 0.05;
double femax = 1e-3;

double xp = 3.5,           //Start location
       w = 0.1,
       h = 0.5,            //Max. height
       hb = 0.;            //Height above bottom floor

#define RATIO (1./772.4) //rho2 = RATIO; two-phase.h
#define MURATIO (17.9e-6/1.12e-3) //mu2 = mu1*MURATIO; two-phase.h
//#define U_X (t > 1./(2.*f_) ? 1e-6 : (1.*s_*pi*f_*sin(2*pi*f_*t)))
//#define U_X (0.)
#define PISTON (w - fabs(x - xp) - (y > h) - (y < hb))
#define g_   9.81 //G.y = -g_; reduced.h
#define SIGMA 0.0728

scalar pstn[];

int main() {
  L0 = 5.45;
  rho1 = 1018.3;
  rho2 = rho1*RATIO;//two-phase.h
  mu1 = 1.12e-3;
  mu2 = mu1*MURATIO;//two-phase.h
  f.sigma = SIGMA;
//  TOLERANCE = 1e-4; 
  G.y = -g_;
//  DT = 1e-5;
  run();
}

u.t[left] = dirichlet(0);
u.t[right] = dirichlet(0);

pstn[left] = dirichlet(0);
pstn[right] = dirichlet(0);

event init (t = 0) {
  if (!restore (file = "restart")) {

  fraction (f, (0.25 - y) - (x > 2.5)); //Water depth $h = 0.25$

  pstn.refine = pstn.prolongation = fraction_refine;
  fraction (pstn, PISTON);

  }
}

event piston (i++; t <= 2.) {
  foreach() 
    foreach_dimension() {
    u.x[] *= (1 - pstn[]);
  }
  boundary ((scalar *){u});
}

event movie (t += 0.01) {
  scalar omega[];
  vorticity (u, omega);

  view (tx = -0.5);
  clear();
  draw_vof ("pstn", fc = {0.2,0.2,0.2});
  draw_vof ("f", filled = 1, fc = {0.1,0.1,0.9});
  begin_mirror ({0,-1});
  cells();
  end_mirror();

  save ("movie.mp4");
}

event cfl (i++) {
    CFL = 0.5;
}

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({f, u}, {pstn}, (double[]){femax, uemax, uemax, uemax}, maxlevel, minlevel, 0);
}
#endif

event end (t = 2.0) {
  fprintf (fout, "i = %d t = %g\n", i, t);
}
