/* This is a lid driven cavity program specifically lid.c on http://basilisk.fr/src/test/lid.c with Re=1000 */ 
// Author: Anvesh 
// This code is primarily written to validate the results with Ghia Et Al.
//The libraries include- 


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "embed.h"



char name[100];
char name_vtk[100];

int maxlevel = 8;  

int main()
{

	//COORDINATES OF LOWER-LEFT CORNER 
	

	origin(-0.5,-0.5);

	// no of grid points 
	
	const face vector muc[] = {1e-3,1e-3};
	mu = muc;

	N=256;


	// maximum timestep
	DT = 0.1;
	// CFL number
	CFL = 0.1;
	
	run(); 



}


	u.t[top] = dirichlet(1); 


	u.t[bottom] = dirichlet(0);
	u.t[left]   = dirichlet(0);
	u.t[right]  = dirichlet(0);








static double energy()
{
  double se = 0.;
 	 if (u.x.face)
    		foreach(reduction(+:se))
      		se += (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.*sq(Delta);
 	 else // centered
    		foreach(reduction(+:se))
      		se += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
 	 return se;
}	


	scalar un[];


event logfile (t += 0.1; i <= 10000) 


{
	  double du = change (u.x, un);
	  if (i > 0 && du < 1e-5)
	    return 1; /* stop */
	  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}







	event outputfile (i += 100) 
{
//	output_matrix (u.x, stdout, N, linear = true);

}


event movies (i += 10 )
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
        output_ppm (u.x, fp, min = -2, max = 2, n = 512);


}




















event profiles (t = end)
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, 0, y));
  fclose (fp);

  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}



















