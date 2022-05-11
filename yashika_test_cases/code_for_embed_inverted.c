#include "embed.h"
#include "navier-stokes/centered.h"
#include "vof.h"
double rho1 = 1., mu1 = 0.01, rho2 = 0.85, mu2 = 0.255, radius = 0.005*0.95;
#include "tension.h"

scalar f[], * interfaces = {f};
face vector alphav[], muv[], av[];
scalar rhov[];

int MAXLEVEL = 10;

u.n[left]  = dirichlet(y<0.005 && y>-0.005 ? 0.0423858 : 0.0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
f[embed]   = 0.;

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
f[bottom]   = 0.;  

int main() {
  size (0.1);

  origin (-0.01, -0.035);

  a = av;
  alpha = alphav;
  rho = rhov;
  mu = muv;
	
  f.sigma = 20.0;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
	
	/*  following 5 lines are for restarting the code in between */
   //~ for (scalar s in {f, u, g})
  //~ s.prolongation = refine_injection;
  //~ restore ("snap-180");
  //~ for (scalar s in {f, u, g})
  //~ s.prolongation = refine_embed_linear;	
	
/* Script for embed function  */	
  vertex scalar phi[];
  foreach_vertex(){ phi[] = ( (y<0.005 && y>-0.005) || ((x>=0.05 && x<=0.08) && y<-0.005) ); }
  boundary ({phi});
  fractions (phi, cs, fs);
  refine (cs[]>0.0 && level<MAXLEVEL);
  
  /* initial condition */
  fraction (f, sqrt(sq(x-0.035)+y*y)<0.95*radius); /* comment this line when you want to restart the simulation*/
  
 //~ boundary(all); /* UNcomment this line when you want to restart the simulation */
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
  foreach_face(y)
    av.y[] -= 800.0;
  boundary ((scalar *){av});
}

event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

event snapshot (t = 0; t += 0.01; t <= 50.) {
  char name1[80];
  sprintf (name1, "snap-%g", t*100);
  p.nodump = false;
  dump(name1);
}

#define IN_REGION ((y<0.005 && y>-0.005) || ((x>=0.05 && x<=0.08) && y<-0.005))
event adapt (i++) {
  scalar f1[];
  scalar region[];
  foreach() 
  f1[] = f[]; 	
  foreach()
  region[] = IN_REGION*noise();
  adapt_wavelet ({region,f}, (double[]){0.01,0.001}, minlevel = 3, maxlevel = MAXLEVEL);
}
