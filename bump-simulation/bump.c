#include"saint-venant.h"


event end (i = 300) 
		{
		
		printf("i =  %d and t=%g  ", i,t )
		
		}

event init( t=0 )
	{
	
	 foreach()
		 h[] = 0.1 + 0.1*exp(-200.*(x*x + y*y ));  


	}	

event graphs(i++) 
	{
	
	stats s = statsf(h);	 
	fprintf( stderr, "%g %g %g \n ", t , s.min , s.max )
	
	
	
	}



event images( i++  )
	{
	
	output_ppm(h); 
	
	}

event adapt (i++) 
	{
	
		adapt_wavelet( {h} , (double []) {4e-3}, maxlevel =8    ); 
	
	
	
	
	}



int main ()
{
 origin(-0.5,-0.5);
init_grid(256); 

run();
return 0; 




}
