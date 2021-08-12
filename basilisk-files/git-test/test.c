// test file that i found on github 
// qcc -O2 test.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm





#define FILTERED //RC
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)) //RC

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

FILE *fp1 ;

#define LEVEL 8 // RC was 4, needs to be bigger to capture the setup

// Dimensional quantities:
#define rhoWater 1000.0
#define rhoOil 917.0
#define muWater 0.001 // approximatley the viscosity of water
#define muOil 0.03 // approximateley the viscosity of oil

#define sig 0.0728  //surface tension of water

#define refLength 0.00001  // height of the domain
#define refVelocity 10.0  // velocity of the top plate

// Dimensionless quantities:
#define rho_ratio (rhoOil/rhoWater)
#define mu_ratio (muOil/muWater)
#define Re (rhoWater*refVelocity*refLength/muWater)  // Reynolds number
#define We (rhoWater*pow(refVelocity,2)*refLength/sig)

FILE *fp_params;

// RC
FILE * fp_stats;

FILE * fp_interface;

u.t[top] = dirichlet(1.0);
u.n[top] = dirichlet(0.0);

u.t[bottom] = dirichlet(0.0);
u.n[bottom] = dirichlet(0.0);

// scalar s1[];

// void draw_frame(char * fname)
// {
// 	cells();
// 	squares("s1");
// 	save("fname.png");
// }

int main() {
  L0 = 8.;
  origin(-L0/2., -0.5);
  periodic(right);
  init_grid (1 << LEVEL);

  // RC
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // I rescale the system so that water has unit density, and the oil properties are just defined in terms of water
  rho1 = 1.0; // water density
  rho2 = rho_ratio; // oil density
  mu1 = 1/Re; // water dynamic viscosity
  mu2 = mu_ratio*mu1; // oil dynamic viscosity
  f.sigma=1/We; // change this later

  DT = 1.0e-2; // RC
  NITERMIN = 1; // default 1
  NITERMAX = 200; // default 100
  TOLERANCE = 1e-4; // default 1e-3

  char params[2000];
  sprintf(params, "params.txt");
  fp_params=fopen(params, "w");
  {
    char params[2000];
	sprintf(params, "params.txt");
	fp_params=fopen(params, "w");
  }

  fprintf(fp_params, "rho1: %g \n", rho1);
  fprintf(fp_params, "rho2: %g \n", rho2);
  fprintf(fp_params, "mu1: %g \n", mu1);
  fprintf(fp_params, "mu2: %g \n", mu2);
  fprintf(fp_params, "sigma: %g \n", f.sigma);
  fprintf(fp_params, "refLength: %g \n", refLength);
  fprintf(fp_params, "refVelocity: %g \n", refVelocity);
  fprintf(fp_params, "Reynolds Number: %g \n", Re);
  fprintf(fp_params, "Weber Number: %g \n", We);
  fclose(fp_params);

  {
    char intName[2000];
    sprintf(intName, "loginterface.dat");
    fp_interface = fopen(intName, "w");
  }
  run();

  fclose(fp_interface);
  fclose(fp_stats);
}

event init(t = 0) {
 
  // initially velocity is 0 everywhere
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }

  // RC cut the top and bottom of domain
  mask (y > 0.5 ? top : none);
  mask (y < -0.5 ? bottom : none);
  
  fraction (f, 0.3-y);

  // boundary conditions
  u.t[top] = dirichlet(1.);
  u.t[bottom] = dirichlet(0.);

}

// event adapt (i++)
// {
//   adapt_wavelet({u, f}, (double[]){1e-2, 1e-2, 1e-2}, LEVEL);
// }
event adapt(i++)
{
	// Order matters for the following (in which the regions called by refine() and unrefine() overlap. If they don't overlap, you can at least switch refine() and unrefine(), idk about other rearrangements.)
	// unrefine(y<0.1);
	// refine(y<0.1&&level<9);
	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, (LEVEL+2), (LEVEL-3));

  // Going fast to get video editing working:
  adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, 7, 4);
  
	// // refine(y>0.4&&x>1.0&&level<=8);
	// unrefine(y<0.1);
	// refine(y<0.1&&level<9);
	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, (LEVEL+2), (LEVEL-3));

	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);
	// unrefine(x>2.0);

	// Order doesn't matter for the following:
	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, 8, 4);
	// refine(y<-0.45 && level<8);

	// Order doesn't matter for the following:
	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);
	// unrefine(y<0.1);

	// Order doesn't matter for the following:
	// The following two lines manually set the refinement level to be higher near the interface (y=0.3):
	// refine(y<0.35 && y>0.25 && level<9);
	// unrefine(y<0.25 && y>0.35);

	// The following line is how I probably want to actually do the AMR for this code:
	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, 8, 4);
	
	// The following line does just u. If you don't put (scalar *) in front of it, an error is thrown since u is a vector but Basilisk expects a scalar
	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);

	// Test Radu suggested:
	// refine(y>0.45&&level<=8);
	// adapt_wavelet((scalar *){f}, (double[]){1e-6}, 8, 4);
	// unrefine(y>0.45&&x>3.0&&level<=8);
}

event end (t = 100) { // RC restricted to 400
  printf ("i = %d t = %g\n", i, t);
}

// RC added this for profiling
/*
event logstats (t += 1.0) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}
*/

event gfsview (t += 1.0) { // RC
    char name_gfs[200];
    sprintf(name_gfs,"Slice-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

// event xmovie (t+=0.1, t<=10)
// {
//  clear();
//  // cells(lc={0.5,0.5,0.5}, lw=0.5);
//  squares("u.x", spread=-1, linear=true, map=cool_warm);
//  draw_vof ("f", lc = {1.0,1.0,1.0}, lw=2);
//  // cells();
//  save ("xmovie.mp4");
// }

// event ymovie (t+=0.1, t<=10)
// {
//  clear();
//  squares("u.y", spread=-1, linear=true, map=cool_warm);
//  draw_vof ("f", lc = {1.0,1.0,1.0}, lw=2);
//  // cells(); // Movie is black when this line is included
//  save ("ymovie.mp4");
// }

// event shearmovie (t+=0.1, t<=10)
// {
//   scalar shear[];
//   foreach()
//   	shear[] = (mu(f[]))*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   // foreach()
//   //   if (y<0.3)
//   //   {
//   //     shear[] = (mu1)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   //   }
//   // foreach()
//   //   if (y>=0.3)
//   //   {
//   //     shear[] = (mu2)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   //   }
//   // foreach()
//   // shear[] = (u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   boundary ({shear});
//   clear();
//   squares("shear", spread=1, linear=true, map=cool_warm);
//   // draw_vof("f", lc = {1.0,1.0,1.0}, lw=2);
//   save("shearMovie.mp4");
//   FILE * fp_shear = fopen("shear", "w");
//   output_field ({shear}, fp_shear, linear=true);
//   fclose(fp_shear);
// }

event interface (t+=1)
{
  char name_interface[1000];
  sprintf(name_interface, "interface_%g.dat", t);

  FILE * fp2 = fopen(name_interface, "w");
  output_facets (f, fp2);
}

event loginterface (t += 1.0) {    
    scalar posX[],posY[];
    position (f, posX, {1, 0});
    position (f, posY, {0,1});
    fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, statsf(f).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max);
    fflush(fp_interface);
}

event profiles (t = 0; t+=1.0; t<=100) // RC restricted the output a little, don't overdo it at first!
{
  FILE * fp = fopen("xprof", "a");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g %g\n", t, y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "a");
  for (double x = -4; x <= 4; x += 0.01)
    fprintf (fp, "%g %g %g\n", t, x, interpolate (u.y, x, 0));
  fclose (fp);
  
  // scalar shear[];
  // foreach()
  //   if (y<0.3)
  //   {
  //     shear[] = mu1.*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // foreach()
  //   if (y>=0.3)
  //   {
  //     shear[] = mu2.*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // // fp=fopen("shearprof", "a");
  // // for (double y = -0.5; y <= 0.5; y += 0.01)
  // //   for (double x = -4; x <= 4; x += 0.01)
  // //     fprintf(fp, "%g ", shear);
  // //   fprintf(fp, "\n");
  // fp=fopen("Xshearprof", "a");
  // for (double y = -0.5; y <= 0.5; y += 0.01)
  //   fprintf (fp, "%g %g %g\n", t, y, shear);
  // fclose (fp);

  // fp=fopen("Yshearprof", "a");
  // for (double x = -4; x <= 4; x += 0.01)
  //   fprintf (fp, "%g %g %g\n", t, x, shear);
  // fclose (fp);
}

event vortmovie (t+=1.0)
{
  scalar omega[];
  vorticity(u, omega);

  view (fov=3, width=2600, height=400);
  clear();
  squares("omega", spread=-1.0, linear=true, map=cool_warm);
  // draw_vof("fs",fc = {0.0,0.0,0.0}, lw=1); 
  // draw_vof("f", lc = {0.0,0.0,0.0}, lw=1); // For some reason Basilisk throws an error if you don't put spaces between 'lc='
  save("Vorticity.mp4");
}
