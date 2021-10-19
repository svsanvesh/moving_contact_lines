/* This is a lid driven cavity program specifically lid.c on http://basilisk.fr/src/test/lid.c with Re=1000 */ 
// Author: Anvesh 
// This code is primarily written to validate the results with Ghia Et Al.
//The libraries include- 


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "embed.h"
#include "view.h"


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



// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(-1.);


u.t[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(5*x);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);











event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);



event visual(i++; t< 3)
{
	draw_vof ("cs", "fs");
	squares ("p");
	box();
	cells(); 
	save ("p.mp4");
	
/*	clear();
	draw_vof ("cs", "fs");
	squares ("u", linear = true, spread = -1);
	box();
	save ("u.mp4");


*/




}
event movies (i += 10 )
{

        sprintf (name_vtk, "data-%g.vtk", t);
        FILE * fpvtk = fopen (name_vtk, "w");
        output_vtk ({u.x,u.y,p},N,fpvtk,1);


}








