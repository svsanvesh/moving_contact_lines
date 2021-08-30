#include "fractions.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "curvature.h"
#include "tension.h"
#include "vtk.h"


#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#define rho(f)  (1./(clamp(f,0,1)*(1./rho1 - 1./rho2) + 1./rho2))

#define rhoL 1000  //density of water
        #define muL 0.001 //viscosity of water
        #define surf 0.072 // surface tension air-water
        #define rhoG 1.225 //density of air
        #define muG 0.0000181 // viscosity of air
#define L 0.0173 //domain size
        #define grav -9.81 // gravitational acceleration
        #define radb 0.000867// radius bubble
        #define ylim 0.002 // 0.00150
        #define LEVEL 8
        #define U0 2.0 // velocity
int MINLEVEL = 5;
int MAXLEVEL = 9;
//scalar f1[],f3[],f2[], *interface = {f1,f3,f2};
  //  origin (L/2., 0.);
// Boundary conditions
        u.n[right] = neumann(0.);
       u.t[left] = dirichlet(0.);
 u.t[right] = dirichlet(0.);
        u.n[left] = neumann(0.);

       char name[100];

//main code starts
       int main() {
        size (L);
  init_grid (1 << LEVEL);
  rho1 = rhoL;
  rho2 = rhoG;
  mu1 = muL;
  mu2 = muG;
  
  refine (sq(x-(L/10) ) + sq(y) <= sq(radb) && level <= LEVEL);
  refine(x<=L/2 && level<=LEVEL);

//origin (L/2., 0.);
  f.sigma = surf;
  DT = 1e-3;
  TOLERANCE = 1e-5;
  run();
}
event init(t=0) {


    //      fraction (f1, sq(x-L/5 ) + sq(y) - sq(radb));
//fraction(f3, 0.015)
//	refine (sq(x-(L/10) ) + sq(y) <= sq(radb) && x<=L/2 && level <= LEVEL);
      foreach(){
        if(sq(x-(L/10) ) + sq(y) < sq(radb))
                f[] = 0;
        else if(x<L/2)
                f[] = 1;
        else
                f[] = 0;
      }
         output_facets(f,stdout);
      boundary(all);
}
/*            event velocity(t=0)
      {
              face vector uv = u;
                      foreach()
                      {
      if(sq(x-(10*radb)) + sq(y)< sq(radb))
              uv.x[] = 0.01;

      else
              uv.x[] = 0;

                      }
 }*/



 
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] = grav;
}
 




event adapt (i++) {
adapt_wavelet ({f,u.x, u.y}, (double[]){1e-4, 5e-3, 5e-3},
                 maxlevel = LEVEL, minlevel = LEVEL - 2);
}
/*event velocity(t +=0.001;t<=0.1){
        double v_x = 0.;
        double vol = 0.;
	
        foreach(){
                v_x +=   u.x[] * dv() ;
                vol +=  dv();
        }
	
        double  V_x = v_x / vol ;
      
	
  FILE * fp_v = fopen ("v.dat", "a");
  fprintf(fp_v,"%.6f %.6f\n",t,V_x);
  fclose(fp_v);
}*/

event logfile (i++; t <= 2.5)
{
  if (i == 0)
    fprintf (ferr, "t \t\t dt \t\t grid->tn\n\n");
  fprintf (ferr, "%.6f \t %.6f \t %ld\n", t, dt, grid->tn);
}
event drop_image (t += 0.0001)
{
  char name[50];
  sprintf (name, "image_%.2f.ppm", t);

  FILE * fp = fopen (name, "w");

  scalar f0[];
//scalar m[];
  foreach()
    f0[] = f[];
  boundary ({f0});

  output_ppm (f0, fp, 512,min = 0,max=1.0, linear = true);
  fclose (fp);
}

event movies (t=t+0.005; t <=0.1)
{
	scalar omega[];
   vorticity (u, omega);
   output_ppm (omega, file = "vort.mp4",512,
              min = -266, max = 266, linear = true);
   sprintf(name, "vort-%g.vtk", t);
   FILE * fp = fopen (name, "w");
   output_vtk({f,u.x,u.y,p},N,fp,1);
   output_ppm (f, file = "f.mp4",512,
                        min = 0, max = 1.0, linear = true);
   sprintf (name, "data-%g.vtk", t);
   FILE * fpvtk = fopen (name, "w");
        output_vtk ({f,u.x,u.y,p},N,fpvtk,1);
}


