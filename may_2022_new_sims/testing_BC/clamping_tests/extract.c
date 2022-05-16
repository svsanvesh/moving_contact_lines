#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "output_vtu_foreach.h" // new file which works for adaptive grid (vtk only works for uniform
#include <sys/stat.h> //  to add directories -- Palas

int main(){
	// Create directories to avoid seg faults -- Palas
	struct stat st = {0};
	char name[80];
	sprintf (name, "data");
	run();
}

event init (t = 0.0){
	char filename[80];
	int pt = 804;

	while(1){
	  sprintf (filename, "dump-%d", pt );

		printf("extracting %s\n",filename);
		fflush(stdout);

		if (!restore (file = filename)){

			printf("files ended");
			exit(0);
		}
		scalar *list;
		bool linear;
		char name2[80];
		sprintf(name2, "snap-%d.vtu",pt);
		FILE *ptr2 = fopen(name2,"w");

		output_vtu_bin_foreach ((scalar *) {f, p}, (vector *) {u}, N, ptr2, false);

	fclose(ptr2);

	pt = pt + 1;
	}
}
