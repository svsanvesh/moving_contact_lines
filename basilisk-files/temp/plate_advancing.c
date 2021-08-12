// This file is used to solve basic flow problem(axisymmetric)


//Created by Amartya 


// Libraries Included
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "utils.h"
#include "two-phase.h"
#include "tension.h"

uf.n[top] = 0.;

scalar f0[];
u.n[left]  = dirichlet(f0[]*(0.00145));
u.t[left]  = dirichlet(0);

p[left]    = neumann(0);
f[left]    = f0[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);


scalar un[];

int maxlevel = 9;   

char name[500];
char name_vtk[500];
char name1[500];
double uemax = 0.1;
int j;


#define Rectangle (min(min(x,0.0007-x),min(0.00047-y,-0.00035+y)))


face vector muc[];
int main() {
  L0 = 0.0035; 
  // origin (-0.5,-L0/2);
  N = 64;  
  rho1 = 997.;
  rho2 = 1.23;
  mu1 = 0.00089;
  mu2 = 0.000018;
  f.sigma = 0.072;
  j = 0;
  run();
}

bid rectangle;
u.t[rectangle] = dirichlet(0);

event init (t = 0) {
  refine (level <= maxlevel*(1. - sqrt(fabs(-Rectangle))));
  fraction (f0, x<0.0001);
  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0});
  vertex scalar phi[];
  foreach_vertex()
      phi[] = -Rectangle;
  boundary ({phi});
  
  if (j == 0) {
    unrefine (y < 0.0035 && level > 5);
    mask(y > 0.000875 ?  top : y > 0.00047 && x < 0.0001 ?  rectangle : Rectangle >  0. ? rectangle : none)

    foreach() {
      f[] = f0[]*(x < 0.0001);
      u.x[] = f[];
    }
    boundary ({f,u.x});
  }
}



// Generating output in VTK format for ParaView


event movies (i+=4;t<=1) {
  
  sprintf (name, "tpf_ns_vor-%g.ppm", t);
  sprintf (name_vtk, "data-%g.vtk", t);
  FILE * fpvtk = fopen (name_vtk, "w");
  FILE * fp = fopen (name, "w");
  output_vtk ({u.x,u.y,p,f},N,fpvtk,1);
  output_ppm (u.x, fp, min = -2, max = 2, n = 512);

  sprintf (name1, "tpf_ns_vf-%g.ppm", t);
  FILE * fp1 = fopen (name1, "w");
  output_ppm (f, fp1, min = -2, max = 2, n = 512);
  

//  scalar omega[]; 
//  vorticity (u,omega);
//  output_ppm (omega, n = 512, file = "movie_ns_5_l9.mp4");
//  output_ppm (f, n = 512, file = "movie_ns_vf_5_l9.mp4");
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, maxlevel);
}

