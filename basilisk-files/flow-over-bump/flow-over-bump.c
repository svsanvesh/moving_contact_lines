//THIS CODE IS USED TO IMPLEMENT A FLOW OVER A BUMP 
//THE CODE IS PRIMAMRY WRITTEN TO TEST THE EMBED.H HEADER AND THE TWO-PHASE.H HEADER.
//Author: ANVESH 
//DATE: 27- 07-2021
//LIBRARIES USED 



#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "view.h"
#include "navier-stokes/perfs.h"

#define H0 0.6
#define QL 1.
#define BA 0.4
#define NU 1e-2
char name[100];
char name_vtk[100];



int maxlevel = 8;              // Maximum mesh refinement
//N=256; 


double hl = H0;

u.n[left] = dirichlet((t < 10. ? t/10. : 1.)*3./2.*QL/hl*(y < hl ? 1. - sq(y/hl - 1.) : 1.));
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);



int main()
{
  size (30);
  N = 32;
  rho1 = 0.001;
  rho2 = 1.;
  mu2 = NU;
  mu1 = mu2*rho1/rho2/10.;
  G.y = - 9.81;
  Z.y = H0;
  run();
}

event init (i = 0)
{
  refine (y < 1.5 && level < 10);
  
  fraction (f, y - H0);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y - BA*exp(- sq(x - 10.)/5.) - 1e-3;
  boundary ({phi});
  fractions (phi, cs, fs);
}


event update_hl (i++)
{
  hl = 0.;
  foreach_boundary (left, reduction(+:hl))
    hl += Delta*f[];
  hl = L0 - hl;
  printf ("%g %g\n", t, hl);
}


event snapshot (i += 10; t <= 700000) {
  p.nodump = false;
  dump();
}

event maxdt (t <= 10.; t += 0.5);

event profiles (t += 5)
{
  fprintf (stderr, "# prof%g\n", t);
  output_facets (f, stderr);
}



/*
// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);
*/



event pictures (t = 700005)
{
  view (fov = 4.04484, tx = -0.498476, ty = -0.0923365, sy = 5,
	bg = {1,1,1},
	width = 1869, height = 390);
  draw_vof ("cs", filled = -1, fc = {1,1,1});
  draw_vof ("f", filled = 1, fc = {1,1,1});
  squares ("u.x", min = -0.5, max = 4, linear = true);
  isoline ("u.x", 0., lc = {1,1,1}, lw = 2);
  save ("u.x.png");

  draw_vof ("cs", filled = -1, fc = {1,1,1});
  draw_vof ("f", filled = 1, fc = {1,1,1});
  squares ("u.y", min = -0.8, max = 0.8, linear = true);
  save ("u.y.png");  
}

// Produce vorticity animation
/*
event movies (i += 10  ; t <=10)
{
        scalar omega[], m[];
        vorticity (u, omega);
        foreach()
                m[] = cs[]  ;
        boundary ({m});
        output_ppm (omega, file = "vort.mp4", box = {{-L0/2., 0.0},{L0/2., 0.1 }},
                        min = -0.5, max = 2.0, linear = true, mask = m);

}
*/

event adapt (i++) {
        adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 5);
}
































