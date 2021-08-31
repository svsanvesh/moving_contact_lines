#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
//#include "navier-stokes/conserving.h"
//#include "adapt_wavelet_leave_interface.h"
#include "embed.h"

// Mesh
#define LEVEL 10
#define lev_init 7
#define neighbor_levels 1


double tEnd = 1;


#define RE    50.
#define ST    25.
#define CA    1e12


//#define mu_water  0.000852
#define mu_air    0.00001846 
#define mu_water 10*mu_air


//#define rho_water 996.59
#define rho_air   1.177
#define rho_water 10.0* rho_air


#define g 9.81


#define k_1 (RE*mu_water)/(2.*rho_water)
#define k_2 (ST*mu_water)/(4.*rho_water*g)

#define L1    2*B
#define LX    10*B
#define LY    10*B
#define B     pow((k_1*k_2),1./3.)


#define u_m pow(k_1,2./3.)/pow(k_2,1./3.)
#define u_p u_m*(3./2.)

 
#define SIGMA 0

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
f[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.t[left] = dirichlet(0.);
u.n[left] = y/B<1. ? dirichlet(-u_p*sq(y)/sq(B)+u_p) : dirichlet(0.);
f[left] = y/B<1. ? dirichlet(1.) : dirichlet(0.);

u.t[embed] = dirichlet(0.);
u.n[embed] = dirichlet(0.);
f[embed] = y/B < 1. ? dirichlet(1.0) : dirichlet(0.);


u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
f[top] = dirichlet(0.);
p[top] = neumann(0.);
pf[top] = neumann(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0.); 
f[bottom] = dirichlet(1.); 
pf[bottom] = neumann(0.);
p[bottom] = neumann(0.);

int main(int argc, char *argv[]) {

size (LX);
origin (0., 0.);
init_grid (1 << lev_init);

rho1 = rho_water, mu1 = mu_water;
rho2 = rho_air;   mu2 = mu_air;
f.sigma = SIGMA;

run();
}

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (LY - y, y); /* Intersection between plane y < LY AND y > 0 */
    if ( x>-1 && x < L1)
    phi[] = difference (phi[], -B +y ); /* Intersection between Lx = 100*B square and B*(LY-B) rectangle */
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  fraction (f, -y + B);

  double iteration = 0;
  do {
  fraction (f, -y + B);
  iteration++;
  } while (adapt_wavelet ({f}, (double[]){1e-7}, LEVEL,lev_init).nf !=0 && 
    iteration<=10);
  
  scalar m[];
  foreach()
  m[] = 1.0;

  foreach() {
  if(x>-1 && x < L1)
  if (y/B>1.0) 
  m[] = -1.;
  }
  boundary({m});
  static FILE * fp = fopen ("InitialState.ppm", "w");
  output_ppm(f, fp, min=0, max=1, mask = m);

  foreach() {

  if(y/B<1) 
  {
    u.x[] = cs[] ? -u_p*sq(y)/sq(B)+u_p: 0.0;
    f[] = cs[] ? 1.0 : 0.0;
  }
  else 
  {
    u.x[] = cs[] ? 0. : 0.; 
    f[] = cs[] ? 0. : 0.;
  }

  }

  boundary ({f,u.x});

}

event adapt(i++){
double uemax = u_p;
double intemax = uemax/20;
 adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){intemax,uemax,uemax}, LEVEL,lev_init, neighbor_levels);
}



event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
  av.x[] += g;
}

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}

event snapshot (t = 0; t<=tEnd; t+=0.001)
{
  char name[80];
  sprintf (name, "dump_%g", t);
  dump (file = name);
}

event plotInterface (t = 0; t<=tEnd; t+=0.001) {

  char name[80];
  sprintf (name, "interface-%f.txt", t);
  FILE* fp = fopen (name,"w");
  output_facets (f, fp);
  }

event output (t = 0; t<=tEnd; t+=0.001) {
  scalar m[];
  foreach()
  m[] = 1.0;

  foreach() {
  if(x>-1 && x < L1)
  if (y/B>1.0) 
  m[] = -1.;
  }
  boundary({m});

  int n1=pow(2,lev_init);
  char name[80];
  sprintf (name, "t_%g.dat",t);
  FILE * fp = fopen (name, "w");
  output_field ({cs, f}, fp,n=n1, linear = true);
  fclose (fp);
  output_ppm (f, file = "f.png", linear = true, box = {{0.0,0.0},{LX,LY}},min=0,max=1,n=n1, mask = m);
  output_ppm (u.x, file = "u.png", linear = true, box = {{0.0,0.0},{LX,LY}},min=0,max=u_p,n=n1, mask = m);
}

event interface (i+= 20) {
  int n1=pow(2,lev_init);
  static FILE * fp1 = fopen ("liquid_sheet.ppm", "w");
  static FILE * fp2 = fopen("grid.ppm","w");
  output_ppm(f, fp1, n=n1, min=0, max=1, box = {{0,0},{LX,LY}});
  scalar l[];
  foreach()
    l[] = level;
  output_ppm(l,fp2, n=n1,min=5, max=LEVEL, box = {{0,0},{LX,LY}});
}

