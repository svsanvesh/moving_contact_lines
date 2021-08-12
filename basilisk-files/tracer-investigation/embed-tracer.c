#include "tracer.h"
#include "embed.h"
#include "navier-stokes/centered.h"

scalar f[];
scalar * tracers = {f};




int main ()
{
	L0 = 8.;
	origin (-0.5, -L0/2.);


	run(); 
}



event init ( t =0 )
{

	foreach_face()
	uf.x[]=0; 



}
/*

event images (i++)
{
	output_ppm (f);
}
*/
