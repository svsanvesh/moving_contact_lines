// The code is to check the rectangular domain created. 
// compile with: qcc -O2 example-channel-flow.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//
//
//
//
//

#include "embed.h"
#include "navier-stokes/centered.h" 
# define DLM_Moving_particle 1
# define NPARTICLES 1
# define ROTATION 1
# define DLM_alpha_coupling 1

//# include "src/dlmfd/DLMFD_reverse_Uzawa.h"b
# include "view.h"
# define LEVEL 5 
# define MAXLEVEL (LEVEL + 4)
#define width 0.5
#define EPS 1e-14

void channel (scalar cs, face vector fs) {

  vertex scalar phi[];

  foreach_vertex()
    phi[] = difference (y - L0/2. + width/2. - EPS, y - L0/2. - width/2. + EPS);
  
  boundary ({phi});
  fractions (phi, cs, fs);
}
# define rhofluid 1. // fluid density
# define grav 0. // gravity acceleration

# define rhosolid 1.1 // solid density
# define diamsolid 0.01 // particle diameter

int U0=10.; 

int main() {

  periodic (right);
  size (1.);
  origin (0., 0.);

  stokes = true;
  /* TOLERANCE = 1e-5; */
  DT = 1e-2;

  N = 1 << LEVEL;
  init_grid (N);
  run ();
}

u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left]   = dirichlet(0.);  // We give no pressure gradient to check for linear profile
pf[left]  = dirichlet(0.0);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow
pf[right]  = dirichlet(0.);


//Dirichlet couditions-


u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(-U0);


u.n[embed] = dirichlet(y > L0/2. ? 1. : 0.);
u.t[embed] = dirichlet(0.);

face vector muv[];

event init (t = 0) {

  foreach() {
    foreach_dimension()
      u.x[] = 0.;
}

channel (cs, fs);

  mu = muv;
  rho = cm;

  }  

event properties (i++) {

  foreach_face()
    muv.x[] = fm.x[];

  boundary ((scalar *){muv});
}


event movie (i=0; i++) {
  
  view (fov = 24, quat = {0,0,0,1}, tx = -L0/2., ty = -L0/2.,
	bg = {1,1,1});
  cells ();
  draw_vof ("cs", "fs", lw = 5);
  squares ("u.x", linear = true, spread = -1, min = 0., max = 1.);
  save ("ux.mp4");
}
 
event profile (t = 5) 
{
  FILE * fp = fopen ("fields", "w");
  foreach()
    fprintf (fp, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  fclose (fp);
}

