//#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"

double rho1 = 1.0, mu1 = 0.01002, rho2 = 0.001204, mu2 = 0.0001825, velocity = 30.0, surface_tension = 72.0;
#include "tension.h"
#include "output_vtu_foreach.h"

int maxlevel = 8;
scalar f[], * interfaces = {f};
face vector alphav[], muv[], av[];
scalar rhov[];

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u.t[right] = dirichlet(0.);
u.n[right] = dirichlet(0.);

u.t[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);

u.t[top] = dirichlet(velocity);
u.n[top] = neumann(0.);

int main(int argc, char* argv[]) {
  size (1.0);

  maxlevel = atoi(argv[1]);
  N = 1 << maxlevel;
	
  origin (0.0,0.0);

  a = av;
  alpha = alphav;
  rho = rhov;
  mu = muv;

  f.sigma = surface_tension;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0){
    //fraction(f,y<1.0);

    //foreach() {	u.x[] = ( y>0.99 ? vel:0.0); }
    //boundary ({u.x});
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

event acceleration (i++) {
face vector av = a;
foreach_face(x)
av.y[] = -981.0;
 boundary ((scalar *){av});
}

event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

event snapshot (t = 0; t += 0.001; t <= 1.0) {
 /* char name1[80];
  sprintf (name1, "snap-%g", t*10);
  p.nodump = false;
  dump (name1);

  char name[80];
  sprintf (name, "snapshot-%g.vtk", t*1000.0);
  FILE *file;
  file  = fopen(name,"w");
  output_vtk ({u,p,f},N,file,1);
  fclose(file);
*/
char name[80];
sprintf (name, "snap-%g.vtu",t*1000);
FILE * fp ;
fp = fopen(name, "w"); 
output_vtu_bin_foreach ((scalar *) {f,p}, (vector *) {u}, N, fp, false);
fclose (fp);    
}
