// This code is used to introduce an interface in the existing couette flow code (couette.c).
// Created by Anvesh 
// Date : 25 th july 2021
//
// libraries used -
#include"navier-stokes/centered.h"
#include"vtk.h"
#include"two-phase.h"



char name[100];
char name_vtk[100];


scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 1., rho2 = 1., mu2 = 1.;
vector h[];

double theta0 = 90;
h.t[bottom] = contact_angle (theta0*pi/180.);

double Reynolds = 20.;
int maxlevel = 7;
face vector muv[];
double H0=0.4, U0=10.;

// Setting the boundary conditions
u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left]    = dirichlet(0.);  // We give no pressure gradient to check for linear profile 
pf[left]   = dirichlet(0.);


u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);  //we give no pressure gradient - couette flow 
pf[right]  = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(U0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(-U0);



event properties (i++)
{
        foreach_face()
         muv.x[] = U0*H0/Reynolds;
}

event init (t = 0)
{

	fraction (f, - (sq(x) + sq(y) - sq(0.5)));
	foreach()
		u.x[] = 0.01;      //(IT DOESNT )  THIS LINE IS COMMENTED TO CHECK IF THE CODE STILL RUNS 




}
int main ()
{
	L0=1.;
	U0 =10;             // Velocity of the bottom plate
        origin (-L0/2, 0.0);  // Origin is at the bottom centre of the box

	run();


}



// This section helps in generating vtk files for postprocessing in paraview. 

event movies (i += 10  ; t <=10)
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









