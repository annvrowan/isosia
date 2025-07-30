#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>


int main(int argc,char** argv) {

  int i,j;

  long maxitt = 1;
  long maxtime = 1;
  long gnr;
  long istep = 0;
  long tstep = 0;
  long ny = 200;
  long nx = 200;

  /*set number of threads*/
  omp_set_num_threads(24);

  printf("This is ptest...");


  /*#pragma omp parallel shared(A,istep,tstep) private(i,j,gnr) firstprivate(maxitt,maxtime,nx,ny)
    {*/

#pragma omp parallel for default(none) private(i,gnr) firstprivate(ny)
  for (i=0;i<ny;i++) {
    gnr = 0.0;
  }
	
  printf("done\n");
  
  exit(1);
  
  
}
