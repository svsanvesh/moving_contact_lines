//#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"

double rho1 = 1.0, mu1 = 0.01002, rho2 = 0.001204, mu2 = 0.0001825;
double radius =  1.0,  a0 = 0.1;
#include "tension.h"
#include "vtk.h"

int maxlevel = 8;
scalar f[], * interfaces = {f};
face vector alphav[], muv[], av[];
scalar rhov[];

/*
     u.t[right] = dirichlet(0.);

u.n[top] = neumann(0.);
p[top]   = dirichlet(0.);
*/

int main(int argc, char* argv[]) {
  size (4.0);

  maxlevel = atoi(argv[1]);
  N = 1 << maxlevel;
	
  origin (-2.0,-2.0);

  a = av;
  alpha = alphav;
  rho = rhov;
  mu = muv;

  f.sigma = 72.0;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0){
  fraction(f,radius - a0*0.5*(5*pow(cos(atan2(y,x)),3)-3*cos(atan2(y,x))) - sqrt(x*x + y*y));
}

#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f)  (clamp(f,0,1)*(mu1 - mu2) + mu2)

event properties (i++) {
  foreach()
    rhov[] = rho(f[])*cm[];
  boundary ({rhov});
  foreach_face () {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  boundary ((scalar *){muv});
}
/*
event acceleration (i++) {
face vector av = a;
foreach_face(x)
av.y[] = -981.0;
 boundary ((scalar *){av});
}
*/

event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

event snapshot (t = 0; t += 0.001; t <= 5.) {
  char name1[80];
  static int pt = 0;
  sprintf (name1, "snap-%d", pt);
  p.nodump = false;
  dump (name1);
	
  char name[80];
  sprintf (name, "snapshot-%d.vtk", pt);
  FILE *file;
  file  = fopen(name,"w");
  output_vtk ({u,p,f},N,file,1);
  fclose(file);	

  pt = pt+1;
}
