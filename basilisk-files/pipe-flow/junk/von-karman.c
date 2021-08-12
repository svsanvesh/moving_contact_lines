


 


#include "embed.h"
#include "navier-stokes/centered.h" 
// #include " navier-stokes/perfs.h" 

#include "tracer.h"



	scalar f[] ;
 	scalar * tracers = {f } ; 
	double Reynolds = 160.  ;
	int maxlevel = 8 ;
	face vector muv[]; 

	 int main () 
 		 { 	
 			L0 =8. ; 
			origin ( -0.5 , -L0/2.0);
			N = 256 ;
			mu = muv ;

			display_control (Reynolds, 10 , 1000);
			
			display_control (maxlevel , 6,  12);
			run();

 
		 }	




	event properties (i++)
	{
		 foreach_face()
    		muv.x[] = fm.x[]*0.125/Reynolds;
	}






	u.n[left]  = dirichlet(1.);
	p[left]    = neumann(0.);
	pf[left]   = neumann(0.);
	f[left]    = dirichlet(y < 0);

	u.n[right] = neumann(0.);
	p[right]   = dirichlet(0.);
	pf[right]  = dirichlet(0.);	



	//u.n[embed] = fabs(y) > 1. ? neumann(0.) : dirichlet(0.);
	//u.t[embed] = fabs(y) > 1. ? neumann(0.) : dirichlet(0.);

	event init (t = 0)
		{
			vertex scalar phi[];
			  foreach_vertex()
			 	{
			    		phi[] = intersection (0.5 - y, 0.5 + y);
			   		phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
		 		}
 	
			 
			 
			 
			 boundary ({phi});
  			fractions (phi, cs, fs);
			foreach()
			u.x[] = cs[] ? 1.0 : 0.; 
		}




		event logfile (i++)
		  	{
       
			fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

			}

		event movies (i += 4; t <= 4.)
		{
			  scalar omega[], m[];
			  vorticity (u, omega);
			  foreach()
			    m[] = cs[] - 0.5;
			  boundary ({m});
			  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
				      min = -10, max = 10, linear = true, mask = m);
			  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
				      linear = false, min = 0, max = 1, mask = m);
		}



		event adapt (i++) 
		{
			  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
		}	


