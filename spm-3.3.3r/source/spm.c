const double pi = 3.14159265358979;
const double s_per_y = 3600.0*24.0*365.25;
const double tiny = 1.0e-16;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

typedef struct {
  double x,y;
  double bed,topice,ice,oldice,newice,snow,phi,isostasy;
  double sedi,dHw,dHs,topsedi,Hgw,Pw,srate,rfail,Tair,Ts,Tb,Ta;
  double maxacc,meltrate,Mb,Ms,accrate,mrate,Kw,tprop,sfac,lee;
  double periglacial_erate,periglacial_erosion,fluvial_erate,fluvial_erosion;
  double hillslope_erate,hillslope_erosion,landslide_erosion;
  double quarrying_rate,quarrying,abrasion_rate,abrasion;
  double dhdx,dbdx,dhdy,dbdy,dtdx,dtdy,alpha2,curv,bcurv,bslope,tslope,hslope;
  double exx,eyy,ezz,exy,hcond,hconv,kappa;
  double te2,te2i,td2,te2r,sxx,syy,sxy,szz,sxz,syz;
  double lb,pb,tn,te,ts,tbx,tby,tbz;
  double vb,vx_b,vy_b,vx_d,vy_d,sliding,deformation,vbres;
  double dsxxdx,dsyydx,dsxydx,dsxxdy,dsyydy,dsxydy;
  double cx[2],cy[2],kc[3];
  double water,Qw,Wc,tau,Qt,streamload;
  int fixflag,ri,rj,di,dj,nparent,streamorder,cnumber,inslide,margin;
  int reciever,nne,ne_i[8],ne_j[8],pa_i[8],pa_j[8];
  double wsurf,lakewater,ndist[8],rdist;
  double Vs[20];


} celltype;

typedef struct {
  double dbdx,dhdx,dtdx,dhdy,alpha2,Kgw;
  double dH,dHw,dHs,kappa;
  double vx,vrx,vx_s,vx_d,vx_b,vz_b,ice,vresx,vx_di,vx_bi;
  double dsxxdx,dsyydx,dsxydx;
  double cx[2],cy[2],kc[3],wx[4];
} hptype;

typedef struct {
  double dbdy,dhdy,dtdy,dhdx,alpha2,Kgw;
  double dH,dHw,dHs,kappa;
  double vy,vry,vy_s,vy_d,vy_b,vz_b,ice,vresy,vy_di,vy_bi;
  double dsxxdy,dsyydy,dsxydy;
  double cx[2],cy[2],kc[3],wy[4];
} vptype;

typedef struct {
  double exy,sxy;
} cornertype;

typedef struct {
  int nx,ny,nc,periodic;
  int nmoni,slidingmode,hydromode;
  int dofluvial,doice,dosliding,doperiglacial,doglacialhydrology,doglacialsedi,dodebrisablation;
  int doglacialerosion,doavalance,dohillslope,dohillslopeerosion,dosediment,coldbased;
  int dolandslide,doisostasy,docelldata,ncelldata,celldata_i[100],celldata_j[100];
  long stressitt,stressitt_count,veloitt,veloitt_count,hydroitt;
  double dx,dy,hmin,ct,L,H,gravity;
  double maxtime,filetime,maxdt,maxdt_ice;
  double meanice,meanele,maxspeed,maxb,maxs,meanice_b;
  double mean_periglacial_erate,mean_dHs_periglacial,mean_dHs_hillslope,mean_dHs_fluvial;
  double mean_quarrying_rate,mean_abrasion_rate,mean_fluvial_erate,mean_hillslope_erate,mean_landslide_erate;
  double meanbed,minbed,maxbed,mean_isostasy,meansedi,meanpw;
} meshtype;

typedef struct {
  double gamma,gamma0,Cs,latentheat,ki,rho,cp;
  double ifac,sfac,vbfac,mf,Ka,ap,qfac;
  double C,L0,minbeta,maxsliding,maxdeformation;
  int maxitt_v,maxitt_s;
} iproptype;

typedef struct {
  int nTemp,nMrate,mtype;
  double lrate,T0,Tsl,maxacc,maxabla,accgrad,ablgrad,avaslope,avacurv,maxslope;
  double qb,sedifrac,Ldebris;
  double **Temp,**Mrate;
} mproptype;

typedef struct {
  double Kgw,a2w,po,tscale;
} hwproptype;

typedef struct {
  int Nc;
  double Ks,sc,Ke,Ls,gamma;
} hproptype;


typedef struct {
  double pr,rho_s,Dg,kw,tau_c,Kt,Ke;
} fproptype;


typedef struct {
  int nHs,nT0;
  double rho_b,rho_s,Ke,Kt;
  double maxsedi,maxice,minslope,minsedi;
  double *Hsv,*T0v;
  double **Ci,**Tr;
} pgproptype;

#include "periglacial.c"
#include "hillslope.c"
#include "hydrology.c"
#include "inout.c"
#include "glacial.c"
#include "fluvial.c"
#include "flexure.c"
#include "thermal.c"

void get_gradients(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh)
{

  int i,j;
  double meanbed = 0.0;
  double meansedi = 0.0;
  double lmaxbed = -10e3;
  double lminbed = 10e3;
  double maxbed = -10e3;
  double minbed = 10e3;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int gmode = 1;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double maxb = (*mesh).maxb;
  double maxs = (*mesh).maxs;


#pragma omp parallel shared(cells,vp,hp,meanbed,minbed,maxbed) private(i,j) firstprivate(nx,ny,dx,dy,maxs,maxb,lminbed,lmaxbed)
  {

    /*update elevations*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].topsedi = cells[i][j].bed + cells[i][j].sedi;
	cells[i][j].topice = cells[i][j].topsedi + cells[i][j].ice;
      }
    }


    /*compute gradients in x dir*/
#pragma omp for schedule(static) nowait
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	/*hp[i][j].dbdx = (cells[i][j].bed-cells[i][j-1].bed)/dx;*/
	hp[i][j].dbdx = (cells[i][j].topsedi-cells[i][j-1].topsedi)/dx;
	hp[i][j].dhdx = (cells[i][j].topice-cells[i][j-1].topice)/dx;
	hp[i][j].dtdx = hp[i][j].dhdx;
	hp[i][j].ice = 0.5*(cells[i][j].ice+cells[i][j-1].ice);
	if (hp[i][j].dbdx < -maxb) hp[i][j].dbdx = -maxb;
	if (hp[i][j].dbdx > maxb) hp[i][j].dbdx = maxb;
	if (hp[i][j].dhdx < -maxs) hp[i][j].dhdx = -maxs;
	if (hp[i][j].dhdx > maxs) hp[i][j].dhdx = maxs;
      }
      hp[i][0].dhdx = hp[i][1].dhdx;
      hp[i][0].dbdx = hp[i][1].dbdx;
      hp[i][0].dtdx = hp[i][1].dtdx;
    }

    /*compute gradients in y dir*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].dbdy = (cells[i][j].topsedi-cells[i-1][j].topsedi)/dy;
	vp[i][j].dhdy = (cells[i][j].topice-cells[i-1][j].topice)/dy;
	vp[i][j].dtdy = vp[i][j].dhdy;
	vp[i][j].ice = 0.5*(cells[i][j].ice+cells[i-1][j].ice);
	if (vp[i][j].dbdy < -maxb) vp[i][j].dbdy = -maxb;
	if (vp[i][j].dbdy > maxb) vp[i][j].dbdy = maxb;
	if (vp[i][j].dhdy < -maxs) vp[i][j].dhdy = -maxs;
	if (vp[i][j].dhdy > maxs) vp[i][j].dhdy = maxs;
      }
    }
    
#pragma omp for schedule(static) nowait
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dhdx = 0.5*(hp[i][j].dhdx+hp[i][j+1].dhdx);
	cells[i][j].dbdx = 0.5*(hp[i][j].dbdx+hp[i][j+1].dbdx);
	cells[i][j].dtdx = 0.5*(hp[i][j].dtdx+hp[i][j+1].dtdx);
      }
    }
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dhdy = 0.5*(vp[i][j].dhdy+vp[i+1][j].dhdy);
	cells[i][j].dbdy = 0.5*(vp[i][j].dbdy+vp[i+1][j].dbdy);
	cells[i][j].dtdy = 0.5*(vp[i][j].dtdy+vp[i+1][j].dtdy);
	cells[i][j].bslope = sqrt(pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0));
	cells[i][j].hslope = sqrt(pow(cells[i][j].dhdx,2.0)+pow(cells[i][j].dhdy,2.0));
	cells[i][j].tslope = sqrt(pow(cells[i][j].dtdx,2.0)+pow(cells[i][j].dtdy,2.0));
	if (gmode > 0) cells[i][j].alpha2 = pow(cells[i][j].dhdx,2.0)+pow(cells[i][j].dhdy,2.0);
	else cells[i][j].alpha2 = 0.0;
      }
    }

#pragma omp for schedule(static)
    for (i=1;i<ny-1;i++) {
      for (j=1;j<nx-1;j++) {
	cells[i][j].curv = (cells[i][j+1].topice-2.0*cells[i][j].topice+cells[i][j-1].topice)/(dx*dx)+(cells[i+1][j].topice-2.0*cells[i][j].topice+cells[i-1][j].topice)/(dy*dy);
	cells[i][j].bcurv = (cells[i][j+1].bed-2.0*cells[i][j].bed+cells[i][j-1].bed)/(dx*dx)+(cells[i+1][j].bed-2.0*cells[i][j].bed+cells[i-1][j].bed)/(dy*dy);
      }
    }


#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	hp[i][j].dhdy = .5*(cells[i][j-1].dhdy+cells[i][j].dhdy);
	if (gmode > 0) hp[i][j].alpha2 = pow(hp[i][j].dhdx,2.0)+pow(hp[i][j].dhdy,2.0);
	else hp[i][j].alpha2 = 0.0;
      }
    }

#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].dhdx = .5*(cells[i-1][j].dhdx+cells[i][j].dhdx);
	if (gmode > 0) vp[i][j].alpha2 = pow(vp[i][j].dhdx,2.0)+pow(vp[i][j].dhdy,2.0);
	else vp[i][j].alpha2 = 0.0;
      }
    }


  /*compute mean, min, and max metrics*/
#pragma omp for schedule(static) reduction(+:meanbed) reduction(+:meansedi)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanbed += cells[i][j].bed;
      meansedi += cells[i][j].sedi;
      if (cells[i][j].bed < lminbed) lminbed = cells[i][j].bed;
      if (cells[i][j].bed > lmaxbed) lmaxbed = cells[i][j].bed; 
    }
  }
#pragma omp critical 
      { 
	if (lmaxbed > maxbed) maxbed = lmaxbed; 
	if (lminbed < minbed) minbed = lminbed; 
      }

  }/*pragma*/

  meanbed /= (double)((*mesh).nc);
  meansedi /= (double)((*mesh).nc);
  (*mesh).meanbed = meanbed;
  (*mesh).minbed = minbed;
  (*mesh).maxbed = maxbed;
  (*mesh).meansedi = meansedi;

}/*get_gradients*/



void initialize(celltype **cells,hptype **hp,vptype **vp,cornertype **cp,meshtype *mesh,hwproptype hwprop,char *path,char file[200],FILE *fp)
{

  int i,j,k;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;


  printf("nx = %d, ny = %d\n",nx,ny);

  /*read bed topography etc*/
  sprintf(file,"%s/input/meshdata.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/meshdata.input\n",path);
    exit(1);
  }
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].bed,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].ice,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].sedi,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].maxacc,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].mrate,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].srate,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].fixflag,sizeof(int),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].phi,sizeof(double),1,fp);
  fclose(fp);

  /*read englacial sediment concentration*/
  if ((*mesh).doglacialsedi > 0) {

    sprintf(file,"%s/input/Vsdata.input",path); /*filename*/
    if ((fp = fopen(file,"rb")) == NULL) {
      printf("Cannot open file : %s/input/Vsdata.input\n",path);
      exit(1);
    }
    for (k=0;k<20;k++) {
      for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].Vs[k],sizeof(double),1,fp);
    }/*k*/
  }/*if*/
    

  /*read cellnumbers for celldata*/
  if ((*mesh).docelldata > 0) {
    sprintf(file,"%s/input/cellnumbers.input",path); /*filename*/
    if ((fp = fopen(file,"rb")) == NULL) {
      printf("Cannot open file : %s/input/meshdata.input\n",path);
      exit(1);
    }
    fread(&(*mesh).ncelldata,sizeof(int),1,fp);
    for (i=0;i<(*mesh).ncelldata;i++) {
      fread(&(*mesh).celldata_i[i],sizeof(int),1,fp);
      fread(&(*mesh).celldata_j[i],sizeof(int),1,fp);
    }
  }

#pragma omp parallel shared(cells,vp,hp,cp,mesh) private(i,j) firstprivate(nx,ny,dx,dy)
  {

    /*initialize geometry and stress*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].te2 = 1.0;
	cells[i][j].exx = 0.0;
	cells[i][j].eyy = 0.0;
	cells[i][j].exy = 0.0;
	cells[i][j].sxx = 0.0;
	cells[i][j].syy = 0.0;
	cells[i][j].sxy = 0.0;
	cells[i][j].x = (double)j*dx + 0.5*dx;
	cells[i][j].y = (double)i*dy + 0.5*dy;
	cells[i][j].topsedi = cells[i][j].bed + cells[i][j].sedi;
	cells[i][j].topice = cells[i][j].topsedi + cells[i][j].ice;
	cells[i][j].oldice = cells[i][j].ice;
	cells[i][j].newice = cells[i][j].ice;
	cells[i][j].snow = 0.0;
	cells[i][j].periglacial_erosion = 0.0;
	cells[i][j].quarrying = 0.0;
	cells[i][j].abrasion = 0.0;
	cells[i][j].fluvial_erosion = 0.0;
	cells[i][j].sliding = 0.0;
	/*cells[i][j].Hgw = cells[i][j].bed;*/
	cells[i][j].Hgw = 0.0;
	cells[i][j].Pw = 0.0;
	cells[i][j].bslope = 0.0;
	cells[i][j].hslope = 0.0;
	cells[i][j].tslope = 0.0;
	cells[i][j].cnumber = -1;
	cells[i][j].isostasy = 0.0;
	cells[i][j].curv = 0.0;
	cells[i][j].Ts = 0.0;
	cells[i][j].Kw = hwprop.Kgw;
	cells[i][j].Mb = 0.0;
	cells[i][j].Ms = 0.0;
	cells[i][j].sfac = 0.0;
	cells[i][j].margin = 0;
	/*for (k=0;k<20;k++) cells[i][j].Vs[k] = 0.0;*/
      }
    }

    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx+1;j++) {
	hp[i][j].vx = 0.0;
	hp[i][j].vx_s = 0.0;
	hp[i][j].vx_d = 0.0;
	hp[i][j].vx_b = 0.0;
	hp[i][j].vz_b = 0.0;
	hp[i][j].dhdx = 0.0;
	hp[i][j].dtdx = 0.0;
	hp[i][j].dbdx = 0.0;
	hp[i][j].dhdy = 0.0;
	hp[i][j].alpha2 = 0.0;
	hp[i][j].dsxxdx = 0.0;
	hp[i][j].dsyydx = 0.0;
	hp[i][j].dsxydx = 0.0;
	hp[i][j].cx[0] = 0.0;
	hp[i][j].cx[1] = 0.0;
	hp[i][j].cy[0] = 0.0;
	hp[i][j].cy[1] = 0.0;
	hp[i][j].kc[0] = 0.0;
	hp[i][j].kc[1] = 0.0;
	hp[i][j].kc[2] = 0.0;
	hp[i][j].wx[0] = 0.0;
	hp[i][j].wx[1] = 0.0;
	hp[i][j].wx[2] = 0.0;
	hp[i][j].wx[3] = 0.0;
	hp[i][j].Kgw = hwprop.tscale*hwprop.Kgw;
	hp[i][j].dH = 0.0;
	hp[i][j].dHw = 0.0;
	hp[i][j].dHs = 0.0;
	hp[i][j].kappa = 0.0;
      }
    }	
#pragma omp for schedule(static)
    for (i=0;i<ny+1;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].vy = 0.0;
	vp[i][j].vy_s = 0.0;
	vp[i][j].vy_d = 0.0;
	vp[i][j].vy_b = 0.0;
	vp[i][j].vz_b = 0.0;
	vp[i][j].dhdy = 0.0;
	vp[i][j].dtdy = 0.0;
	vp[i][j].dbdy = 0.0;
	vp[i][j].dhdx = 0.0;
	vp[i][j].alpha2 = 0.0;
	vp[i][j].dsxxdy = 0.0;
	vp[i][j].dsyydy = 0.0;
	vp[i][j].dsxydy = 0.0;
	vp[i][j].cx[0] = 0.0;
	vp[i][j].cx[1] = 0.0;
	vp[i][j].cy[0] = 0.0;
	vp[i][j].cy[1] = 0.0;
	vp[i][j].kc[0] = 0.0;
	vp[i][j].kc[1] = 0.0;
	vp[i][j].kc[2] = 0.0;
	vp[i][j].wy[0] = 0.0;
	vp[i][j].wy[1] = 0.0;
	vp[i][j].wy[2] = 0.0;
	vp[i][j].wy[3] = 0.0;
	vp[i][j].Kgw = hwprop.tscale*hwprop.Kgw;
	vp[i][j].dH = 0.0;
	vp[i][j].dHw = 0.0;
	vp[i][j].dHs = 0.0;
	vp[i][j].kappa = 0.0;

      }
    }	
    
  }/*pragma*/

  (*mesh).meanice_b = 0.0;
  (*mesh).stressitt_count = 1;
  (*mesh).veloitt_count = 1;
  (*mesh).stressitt = 0;
  (*mesh).veloitt = 0;
  (*mesh).mean_isostasy = 0.0;
  (*mesh).hydroitt = 0;


}/*initialize*/


void get_sealevel_temperature(mproptype *mprop,double time)
{

  /*interpolate sealevel temperature*/
  int i;
  double fac;

  i = 0;
  while ((time > (*mprop).Temp[0][i])&&(i < (*mprop).nTemp-1)) i += 1;
  if (i == 0)
    (*mprop).T0 = (*mprop).Temp[1][0];
  else if (i == (*mprop).nTemp)
    (*mprop).T0 = (*mprop).Temp[1][i-1];
  else {
    fac = (time - (*mprop).Temp[0][i-1])/((*mprop).Temp[0][i]-(*mprop).Temp[0][i-1]);
    (*mprop).T0 = (1.0-fac)*(*mprop).Temp[1][i-1] + fac*(*mprop).Temp[1][i];
  }

}

void get_surface_temperature(meshtype mesh,mproptype mprop,celltype **cells)
{

  int i,j;
  int nx = mesh.nx;
  int ny = mesh.ny;
  double h;

  double lrate = mprop.lrate;
  double T0 = mprop.T0;

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      /*h = cells[i][j].topice - 200.0; if (h < cells[i][j].bed) h = cells[i][j].bed;*/
      h = cells[i][j].topice;
      cells[i][j].Tair = T0 - lrate*h;
      if ((cells[i][j].ice > 10.0)&&(cells[i][j].Tair > 0.0)) cells[i][j].Ts = 0.0;
      else cells[i][j].Ts = cells[i][j].Tair;
  
    }
  }

}



int main(int argc,char** argv) {

  char file[200];
  int i,j,nx,ny,ndtw;
  int nnx[100],nny[100];
  double dx,dy,dt,dtw,maxHgw,maxspeed,tslope,diff,maxdiff;
  FILE *fg;
  FILE *fp = NULL;

  long ncascade,ncatchment;
  long **catchment,**cascade;
  long **slist;
  double **Whelper;

  long *indx;
  double *Tz;
  double **Tcoef;
  double **al;

  celltype **cells;
  hptype **hp;
  vptype **vp;
  cornertype **cp;
  meshtype mesh;
  iproptype iprop;
  mproptype mprop;
  hwproptype hwprop;
  hproptype hprop;
  fproptype fprop;
  pgproptype pgprop;

  double time = 0.0;
  double ftime = 0.0;
  long step = 0;
  long fnr = 0;
  int Nt = 10;


  /*Read path to working directory*/
  char *path = (char*) malloc(500*sizeof(char));
  path = argv[1];

  /*set number of threads*/
  omp_set_num_threads(20);

  printf("******** This is spm-3.2.0.1 **********\n");

  printf("reading input\n");

  /*read input*/
  printf("  mesh\n");
  mesh = read_mesh(path,file,fp);

  /*read iprop*/
  printf("  iprop\n");
  iprop = read_iprop(path,file,fp);

  /*read mprop*/
  printf("  mprop\n");
  mprop = read_mprop(path,file,fp);

  /*read hwprop*/
  printf("  hwprop\n");
  hwprop = read_hwprop(path,file,fp);

  /*read hprop*/
  printf("  hprop\n");
  hprop = read_hprop(path,file,fp);

  /*read fprop*/
  printf("  fproph\n");
  fprop = read_fprop(path,file,fp);

  /*read pgprop*/
  printf("  pgprop\n");
  pgprop = read_pgprop(path,file,fp);

  /*compute and extract*/
  mesh.dx = mesh.L/((double)mesh.nx);
  mesh.dy = mesh.H/((double)mesh.ny);
  nx = mesh.nx;
  ny = mesh.ny;
  dx = mesh.dx;
  dy = mesh.dy;
  dt = mesh.maxdt;
 
  /*allocate memory*/
  cells = (celltype **) malloc(ny*sizeof(celltype*)); for (i=0;i<ny;i++) cells[i] = (celltype*) malloc(nx*sizeof(celltype));
  hp = (hptype **) malloc(ny*sizeof(hptype*)); for (i=0;i<ny;i++) hp[i] = (hptype*) malloc((nx+1)*sizeof(hptype));
  vp = (vptype **) malloc((ny+1)*sizeof(vptype*)); for (i=0;i<ny+1;i++) vp[i] = (vptype*) malloc(nx*sizeof(vptype));
  cp = (cornertype **) malloc((ny+1)*sizeof(cornertype*)); for (i=0;i<ny+1;i++) cp[i] = (cornertype*) malloc((nx+1)*sizeof(cornertype));
  cascade = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) cascade[i] = (long*) malloc(2*sizeof(long));
  catchment = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) catchment[i] = (long*) malloc(2*sizeof(long));
  slist = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) slist[i] = (long*) malloc(2*sizeof(long));
  Whelper = (double **) malloc((ny+4)*sizeof(double*)); for (i=0;i<ny+4;i++) Whelper[i] = (double*) malloc((nx+4)*sizeof(double));

  indx = (long *) malloc(Nt*sizeof(long));
  Tz = (double *) malloc(Nt*sizeof(double));
  Tcoef = (double **) malloc(Nt*sizeof(double*)); for (i=0;i<Nt;i++) Tcoef[i] = (double*) malloc(5*sizeof(double));
  al = (double **) malloc(Nt*sizeof(double*)); for (i=0;i<Nt;i++) al[i] = (double*) malloc(2*sizeof(double));


  printf("Initializing...");
  initialize(cells,hp,vp,cp,&mesh,hwprop,path,file,fp);
  printf("done\n");

  /*surface water drainage patterns*/
  /*initiate_drainage_info(cells,mesh);*/

  printf("Writing output #%d...",0);
  write_output(cells,hp,vp,mesh,path,file,fp,0);
  printf("done\n");

  /*report to status*/
  sprintf(file,"%s/status.dat",path); /*filename*/
  if ((fp = fopen(file,"w")) == NULL) {
    printf("Cannot open file : status.dat.\n");
  }
  else {
    fprintf(fp,"%2.4e %2.2f %ld\n",time,100.0*time/mesh.maxtime,fnr);
    fclose(fp);
  }


  /*create output file for time series*/
  sprintf(file,"%s/output/tseries.dat",path);
  if ((fp = fopen(file,"wt")) == NULL) {
    printf("could not open file for writing\n");
    exit(1);
  }
  fclose(fp);

  compute_leefactor(cells,mesh);

  printf("Enter time loop\n");
    while (time < mesh.maxtime)
      /*          while (step < 1)*/
    {

      /*climate*/
      get_sealevel_temperature(&mprop,time);

      /*ice surface and bed temperatures*/
      get_surface_temperature(mesh,mprop,cells);

      /*compute gradients*/
      get_gradients(cells,hp,vp,&mesh);
      
      /*fluvial transport and erosion*/
      if (mesh.dofluvial > 0) {

	/*update drainage info*/
	if ((step%10 == 0)||(step == 0))  {
	  get_dnetwork(cells,&mesh);
	  get_cascade(cells,mesh,cascade,catchment,&ncascade,&ncatchment); 
	  get_catchment(ncatchment,catchment,cells);
	}

	/*perform fluvial transport and erosion*/
	fluvial_transport(cells,&mesh,fprop,ncascade,cascade,dt);
      }


      /*glacial component*/
      if (mesh.doice > 0) {

	/*if (((step%100 == 0)||(step == 0))&&(mesh.doavalance > 0)) compute_leefactor(cells,mesh);*/

	/*update margin index*/
	get_margin(cells,mesh);

	/*Compute mass balance*/
	ice_accumulation_rate3(cells,mesh,mprop);
	/*ice_accumulation_rate(cells,mesh,mprop);*/
	/*ice_accumulation_rate2(cells,mesh,iprop,mprop,cascade);*/
	mass_balance(cells,&mesh,iprop,mprop,hwprop,dt);

	/*update ice temperatures*/
	if ((step%10 == 0)||(step == 0)) get_thermal(cells,mesh,iprop,mprop,Nt,indx,Tz,Tcoef,al);

	/*glacial hydrology*/
	if (mesh.doglacialhydrology > 0) {

	  if (mesh.hydromode == 1) {
	    glacial_hydrology3(cells,&mesh,hwprop);
	    /*glacial_hydrology4(cells,&mesh,hwprop);*/
	    dtw = dt;
	  }
	  else {
	    maxHgw = 1.0; for (i=0;i<ny;i++) for (j=0;j<nx;j++) if (cells[i][j].Hgw > maxHgw) maxHgw = cells[i][j].Hgw;
	    dtw = 0.3*pow(mesh.hmin,2.0)/(hwprop.tscale*hwprop.Kgw*maxHgw/hwprop.po);
	    if (dtw < dt) ndtw = ceil(dt/dtw); else ndtw = 1; mesh.hydroitt = ndtw; 
	    for (i=0;i<ndtw;i++) glacial_hydrology(cells,hp,vp,&mesh,hwprop,dt/((double)ndtw));
	  }

	}
	else {
	  dtw = dt;
	  for (i=0;i<ny;i++) for (j=0;j<nx;j++) cells[i][j].Pw = 0.6*cells[i][j].ice;
	}
	  

	/*Compute velocities*/
	isosia(cells,hp,vp,cp,&mesh,iprop);
	
      }

	
      if (mesh.doice > 0) {


	/*Compute time step length*/
	dt = mesh.ct*mesh.hmin/(mesh.maxspeed+1.0);
	if (dt > mesh.maxdt) dt = mesh.maxdt;

	/*Update ice thicknesses*/      
      	get_change(cells,hp,vp,&mesh,dt);

	/*perform glacial erosion*/
	if (mesh.doglacialerosion > 0) glacial_erosion(cells,hp,vp,iprop,&mesh,dt);

	/*feed sediment*/
	/*for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    if ((fabs(cells[i][j].x - 35000) < 1000.0)&&(fabs(cells[i][j].y - 14000) < 1000.0)) 
	      cells[i][j].Vs[0] += dt*0.1;
	  }
	  }*/

	if (mesh.doglacialsedi > 0) glacial_sediment_transport(cells,hp,vp,mesh,dt);


      }

      /*hillslope sediment transport and production*/
      if (mesh.dohillslope > 0)
	hillslope(cells,hp,vp,&mesh,hprop,dt);

      /*bedrock landslides*/
      if (mesh.dolandslide > 0) 
	landslide(cells,&mesh,hprop,slist,dt);

      /*periglacial processes*/
      if (mesh.doperiglacial > 0)
	periglacial(cells,hp,vp,&mesh,mprop,pgprop,dt);

      /*flexural isostasy*/
      if ((mesh.doisostasy > 0)&&(step%100 == 0)) 
	isostasy2(cells,&mesh);



      /*ad hoc handling sediment*/
      
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  /*cells[i][j].sedi -= dt*1.0e-3*cells[i][j].sliding;*/
          if (cells[i][j].sedi < 0.0) cells[i][j].sedi = 0.0;
	  /*if (cells[i][j].sedi > 20.0) cells[i][j].sedi = 20.0;*/
	  /*cells[i][j].sedi = 0.0; */
	}
	}



      /*boundary conditions*/      
      for (i=0;i<ny;i++) {
	/*cells[i][0].bed = cells[i][1].bed;
	  cells[i][nx-1].bed = cells[i][nx-2].bed;*/
	cells[i][0].sedi = 0.0;
	/*cells[i][nx-1].sedi = 0.0;*/
	for (j=0;j<nx;j++) {
	  if ((cells[i][j].reciever < 0)||(cells[i][j].fixflag > 0)) {
	    if (cells[i][j].ice > (500.0 - cells[i][j].bed)) cells[i][j].ice = 500.0-cells[i][j].bed;
	    /*if ((cells[i][j].bed + cells[i][j].Pw) > 100.0) cells[i][j].Pw = 100.0 - cells[i][j].bed;
	      if (cells[i][j].Pw < 0.0) cells[i][j].Pw = 0.0;*/
	    /*if (cells[i][j].ice > 0.0) cells[i][j].ice = 0.0;*/
	    if (cells[i][j].ice < 0.0) cells[i][j].ice = 0.0;
	    cells[i][j].sedi = 0.0;
	  }
	  if (cells[i][j].fixflag <= 0) cells[i][j].bed += cells[i][j].srate*dt;
	  /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].sliding > 1.0)) cells[i][j].sedi = 0.0;*/
	}
      }


      step += 1;
      time += dt;
      ftime += dt;
      
      /*printf("time = %4.2f\n",time);*/


      /*output data file*/
      if (ftime >= mesh.filetime) {
	/*if (step%20 == 0) {*/

	fnr += 1;
	/*printf("  Writing output #%ld...",fnr);*/
	write_output(cells,hp,vp,mesh,path,file,fp,fnr);
	/*printf("done\n");*/
	ftime = 0.0;


	/*report to status*/
	sprintf(file,"%s/status.dat",path);
	if ((fp = fopen(file,"w")) == NULL) {
	  printf("Cannot open file : status.dat.\n");
	}
	else {
	  fprintf(fp,"%2.4e %2.2f %ld\n",time,100.0*time/mesh.maxtime,fnr);
	  fclose(fp);
	}

     
      }

      if (step%mesh.nmoni == 0) {

	sprintf(file,"%s/output/tseries.dat",path);
	if ((fg = fopen(file,"at")) == NULL) {
	  printf("could not open file for writing\n");
	  exit(1);
	}
	fprintf(fg,"%4.4e ",time);
	fprintf(fg,"%4.4e ",dt);
	fprintf(fg,"%4.4e ",mesh.meanice);
	fprintf(fg,"%ld ",mesh.veloitt/mesh.veloitt_count);
	fprintf(fg,"%4.4e ",mesh.maxspeed);
	fprintf(fg,"%4.4e ",mesh.mean_periglacial_erate);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_periglacial);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_hillslope);
	fprintf(fg,"%4.4e ",mesh.mean_quarrying_rate);
	fprintf(fg,"%4.4e ",mprop.T0);
	fprintf(fg,"%4.4e ",dtw); 
	fprintf(fg,"%4.4e ",mesh.meanice_b);
	fprintf(fg,"%ld ",mesh.stressitt/mesh.stressitt_count);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_fluvial);
	fprintf(fg,"%4.4e ",mesh.mean_fluvial_erate);
	fprintf(fg,"%4.4e ",mesh.meanbed);
	fprintf(fg,"%4.4e ",mesh.minbed);
	fprintf(fg,"%4.4e ",mesh.maxbed);
	fprintf(fg,"%4.4e ",mesh.mean_hillslope_erate);
	fprintf(fg,"%4.4e ",mesh.mean_landslide_erate);
	fprintf(fg,"%4.4e ",mesh.mean_abrasion_rate);
	fprintf(fg,"%4.4e ",mesh.mean_isostasy);
	fprintf(fg,"%4.4e ",mesh.meansedi);
	fprintf(fg,"%ld ",mesh.hydroitt);
	fprintf(fg,"%4.4e ",mesh.meanpw);
	fprintf(fg,"\n");
	fclose(fg);
	
	mesh.veloitt = 0;
	mesh.veloitt_count = 1;
	mesh.stressitt = 0;
	mesh.stressitt_count = 1;	
      
	/*output celldata*/
	if (mesh.docelldata > 0) {
	  sprintf(file,"%s/output/celldata.dat",path);
	  if ((fg = fopen(file,"at")) == NULL) {
	    printf("could not open file for writing\n");
	    exit(1);
	  }
	  fprintf(fg,"%4.4e ",time);
	  for (i=0;i<mesh.ncelldata;i++) {
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].bed);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].sedi);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].Tair);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].periglacial_erate);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].periglacial_erosion);
	  }
	  fprintf(fg,"\n");
	  fclose(fg);
	}
	
      }
      

    }/*time*/

  for (i=0;i<ny;i++) free(cells[i]); free(cells);
  for (i=0;i<ny;i++) free(hp[i]); free(hp);
  for (i=0;i<ny+1;i++) free(vp[i]); free(vp);
  for (i=0;i<ny+1;i++) free(cp[i]); free(cp);
  for (i=0;i<ny*nx;i++) free(cascade[i]); free(cascade);
  for (i=0;i<ny*nx;i++) free(catchment[i]); free(catchment);
  for (i=0;i<ny*nx;i++) free(slist[i]); free(slist);
  for (i=0;i<ny+4;i++) free(Whelper[i]); free(Whelper);

  exit(1);


}

