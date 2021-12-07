// This is a basilisk code to simulate a elliptic bubble with major axis =5.mm and
// minor axs = 1.5mm and allow it to evolve it with the help of surface tension
//  forces. 
//  Domain size is 15x15(mmxmm)
//
//  Author : Anvesh 
//  date: 6-Dec-2021
//  status: 
//  Libraries used -


#include "navier-stokes/centered.h"
#include "vtk.h"
#include "adapt_wavelet_leave_interface.h"
#include "two-phase.h"
#include "contact.h"
#include "tension.h"


vector h[];  //HEIGHT FUNCTION 
int maxlevel = 9;              // Maximum mesh refinement
char name_vtk[100];             // vtk file name decleration.
double U0;
double theta0 ;


        #define grav  9.81 // gravitational acceleration
        #define rhoL 1000 //density of water
        #define muL 0.001 //viscosity of water
        #define surf  0.072  // surface tension air-water
        #define rhoG 1 //density of air
        #define muG  0.0000181 // viscosity of air
        #define lc 2.7e-3// capillary length 
        #define maj 0.005 // major aaxis
        #define minor 0.0015 // minor axis  
	 



//make sure that the boundary conditions for the face-centered velocity 
//field are consistent with the centered velocity field 
//(this affects the advection term).

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

int main()
{
        L0 = 0.015;            // Size of the square box
        origin (-L0/2, -L0/2);  // Origin is at the bottom centre of the box

                N = 256;
        stokes = true;
//        theta0 = 160*pi/180.0;
        f.sigma = surf;
        f.height = h;
        display_control (maxlevel, 6, 15);

        rho2 = rhoL;   // fluid 2 is given by f = 0.
        mu2 = muL;
        rho1 = rhoG;   // fluid 1 is given by f =1. 
        mu1 = muG;

        run();

}






//Here the code makes sure the refinement of the interface is high. 
event adapt (i++) {

                scalar impose_refine[], f1[];
                foreach(){
                  f1[] = f[];
          }
          boundary({f1});
          adapt_wavelet({f1}, (double[]){1e-2}, maxlevel ,7 );

        }



// Setting the boundary conditions
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);


u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);



u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.0);


// Printing out standard text outputs on the screen
event logfile (i++)
        fprintf (stderr, "%d %g\n", i, t);




char name[80];
// Produce vorticity animation
event movies (i += 100    ; t <= 30)
{
        sprintf (name, "dump_elliptic-%d", i);
        dump (name);
        foreach()
                sprintf (name_vtk, "datax_elliptic-%d.vtk", i);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,mu.x,mu.y,rho,p,f},N,fpvtk,1);
                 scalar omega[];


}

event init (t = 0)
{
       // the interface shape is given here. 
//      fraction (f , -y );

        fraction (f, sq(x/ 0.005 ) + sq(y/0.0015) - sq(1));

        boundary ({f});
}




event videos (i += 1    ; t <= 30)
{

        output_ppm (f, file = "f_elliptic.mp4",8192,
                        min = 0, max = 1.0, linear = true);


}

