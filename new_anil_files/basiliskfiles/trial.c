#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "contact.h"
#include "reduced.h"
#include "tension.h"
#include "output_vtu_foreach.h"

#define MAXLEVEL 7


event logfile2 (i++) {
	if (i == 0)
	fprintf (ferr,"t dt \n");
	fprintf (ferr, "%g %g \n", t, dt);
}

event adapt (i++)
{
    adapt_wavelet({f,u},(double[]){0.01,0.005,0.005,0.005},MAXLEVEL,5);
}

void backup_fields(scalar f,vector u,int nf)
{
    char name[80], subname[80];
    FILE *fp;
    nf > 0 ? sprintf(name,"sol-%4.4d_n%3.3d.vtu",nf,pid()) : sprintf(name,"sol-0_n%3.3d.vtu",pid());
    fp = fopen(name,"w"); 
    output_vtu_ascii_foreach((scalar *){f},(vector *){u},N,fp,false); 
    fclose(fp);

    #if _MPI
        if (pid() == 0)
        {
            nf > 0 ? sprintf(name,"sol-%d.pvtu",nf) : sprintf(name,"sol-0.pvtu");
            nf > 0 ? sprintf(subname,"sol-%4.4d",nf) : sprintf(subname,"sol-0");
            fp = fopen(name,"w"); 
            output_pvtu_ascii((scalar *){f},(vector *){u},N,fp,subname); 
            fclose (fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
}

event dumpfile (t = 0.;t += 0.01;t <= t_final) 
{
    char name[80];
    sprintf(name,"dump-%g",t*100);
    scalar pid[];
    foreach()
    {
        pid[] = fmod(pid()*(npe() + 37),npe());
    }
    boundary({pid});
    dump(name);
}

