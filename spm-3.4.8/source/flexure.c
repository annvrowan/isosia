
#include "./nrutil.h"
#include "./nrutil.c"


void rlft3(double ***data, double **speq, long nn1, long nn2,
	   long nn3, int isign)
{
  void fourn(double data[], long nn[], int ndim, int isign);
  void nrerror(char error_text[]);
  long i1,i2,i3,j1,j2,j3,nn[4],ii3;
  double theta,wi,wpi,wpr,wr,wtemp;
  double c1,c2,h1r,h1i,h2r,h2i;

  if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3)
    nrerror("rlft3: problem with dimensions or contiguity of data array\n");
  c1=0.5;
  c2 = -0.5*isign;
  theta=isign*(6.28318530717959/nn3);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  nn[1]=nn1;
  nn[2]=nn2;
  nn[3]=nn3 >> 1;
  if (isign == 1) {
    fourn(&data[1][1][1]-1,nn,3,isign);
    for (i1=1;i1<=nn1;i1++)
      for (i2=1,j2=0;i2<=nn2;i2++) {
	speq[i1][++j2]=data[i1][i2][1];
	speq[i1][++j2]=data[i1][i2][2];
      }
  }
  for (i1=1;i1<=nn1;i1++) {
    j1=(i1 != 1 ? nn1-i1+2 : 1);
    wr=1.0;
    wi=0.0;
    for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
      for (i2=1;i2<=nn2;i2++) {
	if (i3 == 1) {
	  j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
	  h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
	  h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
	  h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
	  h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
	  data[i1][i2][1]=h1r+h2r;
	  data[i1][i2][2]=h1i+h2i;
	  speq[j1][j2]=h1r-h2r;
	  speq[j1][j2+1]=h2i-h1i;
	} else {
	  j2=(i2 != 1 ? nn2-i2+2 : 1);
	  j3=nn3+3-(i3<<1);
	  h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
	  h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
	  h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
	  h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
	  data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
	  data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
	  data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
	  data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
	}
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  }
  if (isign == -1)
    fourn(&data[1][1][1]-1,nn,3,isign);
}


#define nSWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(double data[], long nn[], int ndim, int isign)
{
  int idim;
  long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  long ibit,k1,k2,n,nprev,nrem,ntot;
  double tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;

  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	  for (i3=i1;i3<=ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    nSWAP(data[i3],data[i3rev]);
	    nSWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
	for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	  for (i2=i1;i2<=ip3;i2+=ifp2) {
	    k1=i2;
	    k2=k1+ifp1;
	    tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
	    tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
#undef nSWAP

void cosfilt(double l3,long nn3,double l2,long nn2,double **qtot,
	     double dflex,double drho,double grav,double **w,double ***load,double **speq) {


  /*solves the 2D elastic flexural equation on a rectangular region using a cosine transform*/
  /*the solution has zero slope on the boundary, perpendicular to the boundary of the domain */

  double fac,f2,f3,f2i2,f3i3,fac1,fac2;
  long i,j,i1=1,i2,i3,nn1=1,isign;

  /*the cosine transform is obtained in the brute force fashion by extendnig the load*/
  /*periodically, enlarging the data set by a factor of 4; waste of time, but what the heck it's fast enough anyway*/

  /*this is moved to global level  
  load=d3tensor(1,nn1,1,2*nn2,1,2*nn3);
  speq=dmatrix(1,nn1,1,2*2*nn2);
  */

  /*copy load into right locations*/
  /*lower left*/
  for (i=1;i<=nn3;i++)
    for (j=1;j<=nn2;j++) load[1][j][i]=qtot[i][j];
  
  /*lower right*/
  for (i=1;i<=nn3;i++)
    for (j=1;j<=nn2;j++) load[1][j][nn3+i]=qtot[nn3-i+1][j];

  /*upper left*/
  for (i=1;i<=nn3;i++)
    for (j=1;j<=nn2;j++) load[1][nn2+j][i]=qtot[i][nn2-j+1];

  /*upper right*/
  for (i=1;i<=nn3;i++)
    for (j=1;j<=nn2;j++) load[1][nn2+j][nn3+i]=qtot[nn3-i+1][nn2-j+1];

  
  /*fft of load*/
  isign=1;
  rlft3(load,speq,nn1,2*nn2,2*nn3,isign);
  
  fac1=nn1*2*nn2*2*nn3/2.0;
  fac2=drho*grav;
  f2=2.0*pi/(2*l2);
  f2*=f2;
  f3=2.0*pi/(2*l3);
  f3*=f3;
   
  for (i2=1;i2<=2*nn2/2;i2++)     /*y*/
    for (i3=1;i3<=2*nn3/2;i3++) {/*x*/
      f2i2=f2*(i2-1)*(i2-1);
      f3i3=f3*(i3-1)*(i3-1);
      fac=dflex*(f2i2*f2i2+f3i3*f3i3+2.0*f2i2*f3i3)+fac2;
      fac=1.0/fac/fac1;
      load[i1][i2][2*i3-1]*=fac;
      load[i1][i2][2*i3  ]*=fac;
    }/*j*/

  for (i2=2*nn2;i2>=2*nn2/2+1;i2--)   /*y*/
    for (i3=1;i3<=2*nn3/2;i3++) {      /*x*/
      f2i2=f2*(2*nn2+1-i2)*(2*nn2+1-i2);
      f3i3=f3*(i3-1)*(i3-1);
      fac=dflex*(f2i2*f2i2+f3i3*f3i3+2.0*f2i2*f3i3)+fac2;
      fac=1.0/fac/fac1;
      load[i1][i2][2*i3-1]*=fac;
      load[i1][i2][2*i3  ]*=fac;
    }/*j*/

  for (i2=1;i2<=2*nn2/2;i2++)       /*y*/
    for (i3=2*nn3/2+1;i3<=2*nn3/2+1;i3++) {/*x*/
      f2i2=f2*(i2-1)*(i2-1);
      f3i3=f3*(i3-1)*(i3-1);
      fac=dflex*(f2i2*f2i2+f3i3*f3i3+2.0*f2i2*f3i3)+fac2;
      fac=1.0/fac/fac1;
      speq[i1][2*i2-1]*=fac;
      speq[i1][2*i2  ]*=fac;
    }/*j*/

  for (i2=2*nn2;i2>=2*nn2/2+1;i2--)   /*y*/
    for (i3=2*nn3/2+1;i3<=2*nn3/2+1;i3++) {      /*x*/
      f2i2=f2*(2*nn2+1-i2)*(2*nn2+1-i2);
      f3i3=f3*(i3-1)*(i3-1);
      fac=dflex*(f2i2*f2i2+f3i3*f3i3+2.0*f2i2*f3i3)+fac2;
      fac=1.0/fac/fac1;
      speq[i1][2*i2-1]*=fac;
      speq[i1][2*i2  ]*=fac;
    }/*j*/


  /*reverse fft of load*/
  isign=-1;
  rlft3(load,speq,nn1,2*nn2,2*nn3,isign);

  /*copy deflection into right locations*/
  for (i=1;i<=nn3;i++)
    for (j=1;j<=nn2;j++) w[i][j]=-load[1][j][i];
  
}/*cosfilt*/

void flexural_isostasy(celltype **cells,meshtype *mesh,iproptype iprop,flexproptype flexprop,double **W)
{

  long i,j; 
  double **loadg;
  double **speq;
  double ***loadt;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double L = (*mesh).L;
  double H = (*mesh).H; 
  double Te = flexprop.Te;
  double rho_i = iprop.rho;
  double rho_r = flexprop.rho_r; 
  double rho_s = flexprop.rho_a;
  double rho_a = flexprop.rho_a;
  double grav = 9.82;
  double Tflex = 1.0e11*pow(Te,3.0)/(12.0*(1.0-0.25*0.25));
  double mean_isostasy = 0.0; 
 
  loadg = dmatrix(1,ny,1,nx);
  speq = dmatrix(1,1,1,2*2*nx);
  loadt = d3tensor(1,1,1,2*nx,1,2*ny);

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      loadg[i+1][j+1] = grav*(cells[i][j].ice*rho_i+cells[i][j].sedi*rho_s-cells[i][j].abrasion*rho_r); 
	/*loadg[i+1][j+1] = grav*(cells[i][j].sedi*rho_s-cells[i][j].abrasion*rho_r); */
      cells[i][j].bed -= cells[i][j].isostasy;
    }
  }

  cosfilt(H,ny,L,nx,loadg,Tflex,rho_a,grav,W,loadt,speq); 

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].isostasy = W[i+1][j+1]; 
      cells[i][j].bed += cells[i][j].isostasy;
      mean_isostasy += W[i+1][j+1];
    }
  }
  mean_isostasy /= (nx*ny);
  (*mesh).mean_isostasy = mean_isostasy; 
 
  free_dmatrix(loadg,1,ny,1,nx);
  free_dmatrix(speq,1,1,1,2*2*nx);
  free_d3tensor(loadt,1,1,1,2*nx,1,2*ny);

}


void flexural_effective_pressure(celltype **cells,meshtype *mesh,iproptype iprop,double **W)
{

  long i,j; 
  double **loadg;
  double **speq;
  double ***loadt;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double L = (*mesh).L;
  double H = (*mesh).H; 
  double Te = 40.0;//thickness of the plate
  double rho_i = iprop.rho;
  double grav = 9.82;
  double Tflex = 1.0e11*pow(Te,3.0)/(12.0*(1.0-0.25*0.25));//flexural rigidity of the plate
  double mean_te = (*mesh).meante; 
 
  loadg = dmatrix(1,ny,1,nx);
  speq = dmatrix(1,1,1,2*2*nx);
  loadt = d3tensor(1,1,1,2*nx,1,2*ny);

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      loadg[i+1][j+1] = rho_i*grav*cells[i][j].te; 
    }
  }

  cosfilt(H,ny,L,nx,loadg,Tflex,rho_i,grav,W,loadt,speq); 

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].te_s = -W[i+1][j+1]; 
    }
  }
 
  free_dmatrix(loadg,1,ny,1,nx);
  free_dmatrix(speq,1,1,1,2*2*nx);
  free_d3tensor(loadt,1,1,1,2*nx,1,2*ny);

}


void isostasy2(celltype **cells,meshtype *mesh)
{

  int i,j;
  double load,iso,mean_isostasy;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double g = 10.0;
  double rho_a = 3200.0;
  double rho_r = 2800.0;
  double rho_s = 2500.0;
  double rho_i = 950.0;

  mean_isostasy = 0.0;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      
      /*load = -rho_i*g*cells[i][j].ice+rho_r*g*(cells[i][j].periglacial_erosion+cells[i][j].abrasion+cells[i][j].quarrying+cells[i][j].fluvial_erosion+cells[i][j].landslide_erosion+cells[i][j].hillslope_erosion);*/
      load = -rho_i*(cells[i][j].ice-cells[i][j].iniice)+rho_r*(cells[i][j].periglacial_erosion+cells[i][j].abrasion+cells[i][j].quarrying+cells[i][j].hillslope_erosion);
      iso = load/rho_a;				   
      mean_isostasy += iso;

    }
  }
  mean_isostasy /= ny*nx;

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      cells[i][j].bed -= cells[i][j].isostasy;      
      cells[i][j].isostasy = mean_isostasy;
      cells[i][j].bed += cells[i][j].isostasy;      


    }
  }

  (*mesh).mean_isostasy = mean_isostasy;


}

void isostasy3(celltype **cells,meshtype *mesh,iproptype iprop,flexproptype flexprop,double *W0)
{


  int i,j,k;
  double qm,wf,dd;
  
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx; 
  double L = (*mesh).L; 
  double Te = flexprop.Te;
  double rho_i = iprop.rho;
  double rho_r = flexprop.rho_r; 
  double rho_s = flexprop.rho_s;
  double rho_a = flexprop.rho_a;
  double grav = 9.82;
  double D = 1.0e11*pow(Te,3.0)/(12.0*(1.0-0.25*0.25));
  double alpha = pow(4.0*D/(rho_a*grav),0.25);
  double wfm = 0.0; 
 
#pragma omp parallel shared(cells,mesh,wfm,W0) private(i,j,k,qm,wf,dd) firstprivate(nx,ny,dx,L,alpha,rho_i,rho_r,D)
  {

  /*loop grid to find average load in columns*/
#pragma omp for schedule(static)
  for (j=0;j<nx;j++) {
    qm = 0.0;
    for (i=0;i<ny;i++) {
      qm += -rho_i*(cells[i][j].ice-cells[i][j].iniice)+rho_r*(cells[i][j].periglacial_erosion+cells[i][j].abrasion+cells[i][j].quarrying+cells[i][j].hillslope_erosion);
    }
    qm *= grav*dx/((double)ny);
    W0[j] = qm*pow(alpha,3.0)/(8.0*D);
  }

  /*loop 1D grid to sum contributions*/

#pragma omp for schedule(static) reduction(+:wfm)
  for (j=0;j<nx;j++) {
    wf = 0.0;
    for (k=0;k<nx;k++) {
      dd = fabs(cells[0][k].x-cells[0][j].x)/alpha;
      wf += W0[k]*exp(-dd)*(cos(dd)+sin(dd));
      dd = fabs(2.0*L - cells[0][k].x-cells[0][j].x)/alpha; /*mirror at right edge*/
      wf += W0[k]*exp(-dd)*(cos(dd)+sin(dd));
    }
    /*add to grid*/
    for (i=0;i<ny;i++) {
      cells[i][j].bed -= cells[i][j].isostasy;      
      cells[i][j].isostasy = wf;
      cells[i][j].bed += cells[i][j].isostasy;      
    }
    wfm += wf;
  }

  }/*pragma*/


/*save mean isostasy to mesh*/    
  (*mesh).mean_isostasy = wfm/((double)nx);

}

