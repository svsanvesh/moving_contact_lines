#include "navier-stokes/centered.h"
#include "embed.h"

FILE *fp1 ;

#define LEVEL 5
#define Re 100

u.t[top] = dirichlet( 1.0 );
u.n[top] = dirichlet( 0.0 );

u.t[bottom] = dirichlet(0.0);
u.n[bottom] = dirichlet( 0.0 );

int main() {
  L0 = 8.;
  origin(-L0/2., -0.5);
  init_grid (1 << LEVEL);
  run();
}

event init(t = 0) {
  periodic(right);
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection(0.5-y, 0.5+y);
  }
  boundary({phi});
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

event end (t = 400) {
  printf ("i = %d t = %g\n", i, t);
}

event profiles (t = end)
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}
