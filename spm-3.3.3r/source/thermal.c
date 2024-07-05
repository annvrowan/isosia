#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20


void tridag(double *a,double *b,double *c,double *r,double *gam,double *u,long N)
{

  long j;
  double bet;
  double mytiny = 1.0e-36;

  if (b[0] == 0.0) printf("Error 1 in tridag");
  u[0] = r[0]/(bet = b[0]);
  for (j=1;j<N;j++) {
    gam[j] = c[j-1]/(bet+mytiny);
    /*    gam[j] = c[j-1]/(bet);*/
    bet = b[j] - a[j]*gam[j];
    if (bet == 0.0) printf("Error 2 in tridag\n");
    u[j] = (r[j]-a[j]*u[j-1])/(bet+mytiny);
      /*u[j] = (r[j]-a[j]*u[j-1])/(bet);*/
  }
  for (j=(N-2);j>=0;j--) u[j] -= gam[j+1]*u[j+1];

}

void bandec(double **a,int n,int m1,int m2,double **al,long *indx,double *d)
{
  unsigned int i,j,k,l;
  int mm;
  double dum;

  mm = m1+m2+1;
  l = m1;
  for (i=1;i<=m1;i++) {
    for (j=m1+2-i;j<=mm;j++) a[i-1][j-l-1] = a[i-1][j-1];
    l--;
    for (j=mm-l;j<=mm;j++) a[i-1][j-1] = 0.0;
  }
  *d = 1.0;
  l = m1;
  for (k=1;k<=n;k++) {
    dum = a[k-1][0];
    i = k;
    if (l < n) l++;
    for (j=k+1;j<=l;j++) {
      if (fabs(a[j-1][0]) > fabs(dum)) {
	dum = a[j-1][0];
	i = j;
      }
    }
    indx[k-1] = i;
    if (dum == 0.0) a[k-1][0] = TINY;
    if (i != k) {
      *d = -(*d);
      for (j=1;j<=mm;j++) SWAP(a[k-1][j-1],a[i-1][j-1]);
    }
    for (i=k+1;i<=l;i++) {
      dum = a[i-1][0]/a[k-1][0];
      al[k-1][i-k-1] = dum;
      for (j=2;j<=mm;j++) a[i-1][j-2] = a[i-1][j-1]-dum*a[k-1][j-1];
      a[i-1][mm-1] = 0.0;
    }
  }
}

void banbks(double **a,int n,int m1,int m2,double **al,long *indx,double *b)
{
  unsigned long i,k,l;
  int mm;
  double dum;

  mm = m1+m2+1;
  l = m1;
  for (k=1;k<=n;k++) {
    i = indx[k-1];
    if (i != k) SWAP(b[k-1],b[i-1]);
    if (l < n) l++;
    for (i=k+1;i<=l;i++) b[i-1] -= al[k-1][i-k-1]*b[k-1];
  }
  l = 1;
  for (i=n;i>=1;i--) {
    dum = b[i-1];
    for (k=2;k<=l;k++) dum -= a[i-1][k-1]*b[k+i-2];
    b[i-1] = dum/a[i-1][0];
    if (l < mm) l++;
  }
}

void findT(long N,double H,double Ts,double Tm,double qb,double ki,double rhoi,double cpi,double Li,double hconv,double hcond,double sliding,double sbed,double **Tcoef,double **al,long *indx,double *Tz,double *Tb,double *Ta,double Ms,double *Mb,double *sfac)
{

  long i;
  double dz = H/((double)(N-1));
  double qeff,qm,vz,fac,Tb_trial,d;
  double vfac = rhoi*cpi/s_per_y;
  double dT = 3.0;

  /*Ms is surface mass balance, Ms > 0 is accumulation, M < 0 is ablation*/
  /*Mb is basal melting, Mb < 0*/
  (*Mb) = 0.0;

  /*vertical velocity gradient*/ 
  double dvzdz = -(Ms + (*Mb))/H;

  if (Ms < 0.0) Ms = 0.0;

  /*assemble tridiagonal system*/
  /*z axis pointing upwards*/
  for (i=1;i<(N-1);i++) {

    /*vertical velocity*/
    fac = (double)((double)i/((double)(N-1)));
    vz = -fac*Ms + (1.0-fac)*(*Mb);

    /*printf("fac = %4.6f\n",fac);*/

    Tcoef[i][0] = ki/(dz*dz) + .5*vfac*vz/dz; /*prediagonal*/
    Tcoef[i][1] = -2.0*ki/(dz*dz) - 0.0*vfac*dvzdz; /*diagonal*/
    Tcoef[i][2] = ki/(dz*dz) - .5*vfac*vz/dz; /*postdiagonal*/
    Tcoef[i][3] = hconv + hcond; /*load vector*/

  }
  
  /*fix basal heat flow*/
  Tcoef[0][0] = 0.0;
  Tcoef[0][1] = -1.0;
  Tcoef[0][2] = 1.0;
  Tcoef[0][3] = -qb*dz/ki;

  /*fix surface temperature*/
  Tcoef[N-1][0] = 0.0;
  Tcoef[N-1][1] = 1.0;
  Tcoef[N-1][2] = 0.0;
  Tcoef[N-1][3] = Ts;

  /*solve the tridiagonal system*/
  /*tridag(Tcoef[0],Tcoef[1],Tcoef[2],Tcoef[3],Tcoef[4],Tz,N);*/
  for (i=0;i<N;i++) Tz[i] = Tcoef[i][3];
  bandec(Tcoef,N,1,1,al,indx,&d);
  banbks(Tcoef,N,1,1,al,indx,Tz);


  /*estimated basal temperature*/
  Tb_trial = Tz[0]; 

  /*cold based*/
  if (Tb_trial < Tm) {

    /*basal temperature*/
    (*Tb) = Tb_trial;

    /*melting rate*/
    (*Mb) = 0.0;

    /*sliding switch*/
    (*sfac) = 0.0;
    
  }

  /*warm based*/
  else {

    /*assemble tridiagonal system*/
    /*z axis pointing upwards*/
    for (i=1;i<N-1;i++) {
      
      /*vertical velocity*/
      fac = (double)((double)i/((double)(N-1)));

      vz = -fac*Ms + (1.0-fac)*(*Mb);
      
      Tcoef[i][0] = ki/(dz*dz) + .5*vfac*vz/dz; /*prediagonal*/
      Tcoef[i][1] = -2.0*ki/(dz*dz) - 0.0*vfac*dvzdz; /*diagonal*/
      Tcoef[i][2] = ki/(dz*dz) - .5*vfac*vz/dz; /*postdiagonal*/
      Tcoef[i][3] = hconv + hcond; /*load vector*/
      
    }

    /*fix basal temperature*/
    Tcoef[0][0] = 0.0;
    Tcoef[0][1] = 1.0;
    Tcoef[0][2] = 0.0;
    Tcoef[0][3] = Tm;

    /*fix surface temperature*/
    Tcoef[N-1][0] = 0.0;
    Tcoef[N-1][1] = 1.0;
    Tcoef[N-1][2] = 0.0;
    Tcoef[N-1][3] = Ts;

    /*solve the tridiagonal system*/
    /*tridag(Tcoef[0],Tcoef[1],Tcoef[2],Tcoef[3],Tcoef[4],Tz,N);*/
    for (i=0;i<N;i++) Tz[i] = Tcoef[i][3];
    bandec(Tcoef,N,1,1,al,indx,&d);
    banbks(Tcoef,N,1,1,al,indx,Tz);
    
    /*Basal temperature (should be Tm)*/
    (*Tb) = Tz[0]; 

    /*heat for melting*/
    qeff = ki*(Tz[0]-Tz[1])/dz;
    qm = qb - qeff; 

    /*melting rate*/
    (*Mb) = -s_per_y*qm/(rhoi*Li);

    /*frictional heating*/
    (*Mb) -= sliding*fabs(sbed)*9.81/Li;

    if ((*Mb) > 0.0) (*Mb) = 0.0;

    /*sliding switch - average heat flow*/
    /*(*sfac) = qm/qb;*/
    /*if ((*sfac) < 0.0) (*sfac) = 0.0;
      if ((*sfac) > 1.0) (*sfac) = 1.0;*/




    /*!!!!!!!!!!*/
    /*(*sfac) = 1.0;*/


  }


  /*summer access of water*/
  /*  if ((Ts > -5.0)&&(H < 500.0)) (*sfac) += (Ts + 5.0)*(500.0 - H)/(500.0);
  if ((*sfac) < 0.0) (*sfac) = 0.0;
  if ((*sfac) > 1.0) (*sfac) = 1.0;*/

  /*  if ((*Tb) > 0.0) (*sfac) = 1.0;
  else if ((*Tb) > -1.0) (*sfac) = (*Tb) + 1.0;
  else (*sfac) = 0.0; 
  */

  /*(*sfac) = 0.5*(erf(((*Tb)-(Tm-2.0*dT))/dT)+1.0);*/
  (*sfac) = 0.5*(erf((Tb_trial-Tm)/dT)+1.0);
  if ((*sfac) < 0.0) (*sfac) = 0.0;
  if ((*sfac) > 1.0) (*sfac) = 1.0;

  
  /*average temperature*/
  (*Ta) = .5*(Tz[0]+Tz[N-1]);
  for (i=1;i<(N-1);i++) (*Ta) += Tz[i];
  (*Ta) /= (double)(N-1);

}/*findT*/

void get_thermal(celltype **cells,meshtype mesh,iproptype iprop,mproptype mprop,int N,long *indx,double *Tz,double **Tcoef,double **al)
{

  int i,j;
  double ice,Ts,Tm,hcond,hconv,sliding,ts,Ms,Mb,sfac;
  double Ta,Tb;

  int nx = mesh.nx;
  int ny = mesh.ny;
  double dx = mesh.dx;
  double dy = mesh.dy;

  double ki = iprop.ki;
  double rhoi = iprop.rho;
  double cpi = iprop.cp;
  double Li = iprop.latentheat;
  double qb = mprop.qb;
 
  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      ice = cells[i][j].ice; /*ice thickness*/
      Ts = cells[i][j].Ts; if (Ts > 0.0) Ts = 0.0; /*ice surface temperature*/

      if (ice > 10.0) {

	Tm = -8.7e-4*ice; /*basal melting point*/
	
	hcond = 0.0;
	hconv = 0.0;

	/*
	if ((i>0)&&(i<ny-1)&&(j>0)&&(j<nx-1)) {
	  hcond = -ki*((cells[i][j-1].Ts-2.0*cells[i][j].Ts+cells[i][j+1].Ts)/(dx*dx)+(cells[i-1][j].Ts-2.0*cells[i][j].Ts+cells[i+1][j].Ts)/(dy*dy));
	  hconv = rhoi*cpi/s_per_y*((cells[i][j].vx_d+cells[i][j].vx_b)*(cells[i][j+1].Ts-cells[i][j-1].Ts)/(2.0*dx)+(cells[i][j].vy_d+cells[i][j].vy_b)*(cells[i+1][j].Ts-cells[i-1][j].Ts)/(2.0*dy));

	}
	else {
	  hcond = 0.0;
	  hconv = 0.0;
	}
	*/
	sliding = cells[i][j].sliding;
	ts = cells[i][j].ts;
	Ms = cells[i][j].Ms;
	Mb = cells[i][j].Mb;
	
	/*compute steady state thermal profile*/
	findT(N,ice,Ts,Tm,qb,ki,rhoi,cpi,Li,0.0,0.0,sliding,ts,Tcoef,al,indx,Tz,&Tb,&Ta,Ms,&Mb,&sfac);
	
	/*save properties to grid*/
	cells[i][j].Mb = Mb;
	cells[i][j].Tb = Tb;
	cells[i][j].Ta = Ta;
	cells[i][j].sfac = sfac;
	cells[i][j].hconv = hconv;
	cells[i][j].hcond = hcond;
	
	if (fabs(Tb) > 100.0) {
	  printf("  Ms = %2.2e, Mb = %2.2e, Ts = %2.2e, Tm = %2.2e, Tb = %2.2e, ice = %2.2e, sliding = %2.2e, ts = %2.2e\n",Ms,Mb,Ts,Tm,Tb,ice,sliding,ts);
	}


      }
      else {

	/*save properties to grid*/
	cells[i][j].Mb = 0.0;
	cells[i][j].Tb = Ts;
	cells[i][j].Ta = Ts;
	cells[i][j].sfac = 0.0;
	cells[i][j].hconv = 0.0;
	cells[i][j].hcond = 0.0;

      }

    }
  }

}

