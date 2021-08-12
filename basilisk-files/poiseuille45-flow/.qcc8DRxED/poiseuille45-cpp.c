@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
@include "_boundarydecl.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "poiseuille45.c"
//This code is wriiten to simulate a poiseuille flow in an incline position. 
//Author: Anvesh 
//
//
//
//
////We need to link various libraries for *******view.h******** header.
//
//The method to compile the code is --   qcc -Wall -O2 [program].c  -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//

//Libaries used:
//
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "vtk.h"

char name[100];
char name_vtk[100];

int main()
{
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);

  stokes = true;
  TOLERANCE = 1e-7;
N=64; 
//  for (N = 16; N <= 64; N *= 2)
    run();
}



scalar un[];


#define WIDTH 0.5
#define EPS 1e-14

event init(t=0)
{

const face vector g[] = {1.,1.};
  a = g;
  mu = fm;

   vertex scalar phi[];
  foreach_vertex()
    phi[] = difference (union (y - x - EPS, x - y - 0.5 + EPS),	y - x - 0.5 + EPS);
  boundary ({phi});
  fractions (phi, cs, fs);

   u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
    
  for (scalar s in {u})
    s.third = true;

    foreach()
    un[] = u.y[];

}


event logfile (t += 0.01; i <= 1000) 
{
  double du = change (u.y, un);
  if (i > 0 && du < 1e-7)
    return 1; /* stop */
}

event profile (t = end) {
  printf ("\n");
  foreach()
    fprintf (stdout, "%g %g %g \n",i,  x, y);
    //fprintf (stdout, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  scalar e[];
  foreach() 
  {
    double x1 = y - x;
    e[] = u.x[] - 0.25*(sq(WIDTH/2.) - sq(x1 >= 0. ? x1 - 0.25 : x1 + 0.75));
  }
  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
  draw_vof ("cs", "fs");
  squares ("u.x", linear = true, spread = -1);
  save ("u.x.png");
}





event movies (i += 10  ; t <=1000)
{
        scalar omega[], m[];
        vorticity (u, omega);
        foreach()
                m[] = cs[]  ;
        boundary ({m});
        output_ppm (omega, file = "vort.mp4", box = {{-0.5, -0.5},{0.5, 0.5 }},
                        min = -0.5, max = 2.0, linear = true, mask = m);

        sprintf (name, "vort-%g.ppm", t);
        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        FILE * fp = fopen (name, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);
        output_ppm (u.x, fp,min= -10 , max =10, n = 512);


}









#endif
