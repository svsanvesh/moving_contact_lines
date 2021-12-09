// This is a basilisk code to simulate a elliptic bubble with major axis =5.mm and
// minor axs = 1.5mm and allow it to evolve it with the help of surface tension
//  forces. 
//  Domain size is 15x15(mmxmm)
//
//  Author : Anvesh 
//  date: 8-Dec-2021
//  //FOR COMPILING : qcc   elliptic_bubble_no_stokes.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

//  status: 
//  Libraries used -

#include "view.h"
#include "navier-stokes/centered.h"
#include "vtk.h"
#include "adapt_wavelet_leave_interface.h"
#include "two-phase.h"
#include "tension.h"


vector h[];  //HEIGHT FUNCTION 
int maxlevel = 9;              // Maximum mesh refinement
char name_vtk[100];             // vtk file name decleration.


        #define rhoL 1000 //density of water
        #define muL 0.001 //viscosity of water
        #define surf  0.072  // surface tension air-water
        #define rhoG 1 //density of air
        #define muG  0.0000181 // viscosity of air
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

        N = 1 << 5;
        stokes = true;
        f.sigma = surf;
        display_control (maxlevel, 6, 15);

        rho1 = rhoL;   // fluid 1 is given by f =1. 
        mu1 = muL;
        rho2 = rhoG;   // fluid 2 is given by f = 0.
        mu2 = muG;

        run();

}






//Here the code makes sure the refinement of the interface is high. 
event adapt (i += 5) {
  adapt_wavelet ({f}, (double[]){1e-3}, maxlevel );
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
event movies (i += 1000    ; t <= 5 )
{
        sprintf (name, "dump_elliptic_stokes-%d", i);
        dump (name);
        foreach()
                sprintf (name_vtk, "datax_elliptic_stokes-%g.vtk", t);
                FILE * fpvtk = fopen (name_vtk, "w");
                output_vtk ({u.x,u.y,mu.x,mu.y,rho,p,f},N,fpvtk,1);


}

event init (t = 0)
{
       // the interface shape is given here. 

        fraction (f, sq(x/ 0.005 ) + sq(y/0.0015) - sq(1));

        boundary ({f});
}




event videos ( t+=0.00001    ; t <=5 )
{

        output_ppm (f, file = "f_elliptic_stokes.mp4",1024,
                        min = 0, max = 1.0, linear = true);
        clear();
        draw_vof ("f");
        cells();
        box();
        save ("fscalar.mp4");



/*
        clear();
        vectors ("u", scale = .0000025, lc = {0, 1, 0}, lw = .8);
        box ();
        save ("u.mp4");
*/


}

