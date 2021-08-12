// This code is used to introduce an interface in the existing couette flow code (couette.c).
// Created by Anvesh 
// Date : 25 th july 2021
//
// libraries used -
//#include"embed.h"
#include"navier-stokes/centered.h"
#include"vtk.h"
//#include"vof.h"
//#include"contact.h"
//#include"tension.h"
#include"fractions.h"



char name[100];
char name_vtk[100];


scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 1., rho2 = 1., mu2 = 1.;

int main ()
{


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









