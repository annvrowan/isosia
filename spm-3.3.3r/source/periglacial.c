

double get_periglacial_erosion_rate(double Hs,double T0,double bslope,pgproptype pgprop)
{

  /*uses bilinear interpolation*/

  int i,j;
  double erate=0.0,ifac,jfac;

  int nHs = pgprop.nHs;
  int nT0 = pgprop.nT0;
  double minslope = pgprop.minslope;
  double minsedi = pgprop.minsedi;

  if ((bslope < minslope)&&(Hs < minsedi)) Hs = minsedi; 
  if (bslope > 0.7) Hs = 0.0;

  i=0; j=0;
  while ((pgprop.Hsv[i] < Hs)&&(i < nHs-1)) i += 1;
  while ((pgprop.T0v[j] < T0)&&(j < nT0-1)) j += 1;

  if ((i == 0)||(i == nHs-1)) {
    if ((j == 0)||(j == nT0-1)) erate = pgprop.Ci[i][j];
    else { 
      jfac = (T0-pgprop.T0v[j-1])/(pgprop.T0v[j]-pgprop.T0v[j-1]);
      erate = (1.0-jfac)*pgprop.Ci[i][j-1] + jfac*pgprop.Ci[i][j];
    }
  }
  else {
    if ((j == 0)||(j == nT0-1)) {
      ifac = (Hs-pgprop.Hsv[i-1])/(pgprop.Hsv[i]-pgprop.Hsv[i-1]);
      erate = (1-ifac)*pgprop.Ci[i-1][j]+ifac*pgprop.Ci[i][j];
    }
    else {
      ifac = (Hs-pgprop.Hsv[i-1])/(pgprop.Hsv[i]-pgprop.Hsv[i-1]);
      jfac = (T0-pgprop.T0v[j-1])/(pgprop.T0v[j]-pgprop.T0v[j-1]); 
      erate = pgprop.Ci[i-1][j-1]*(1.0-ifac)*(1.0-jfac)+pgprop.Ci[i-1][j]*jfac*(1.0-ifac)+pgprop.Ci[i][j-1]*ifac*(1-jfac)+pgprop.Ci[i][j]*ifac*jfac;
    }
  }

  erate *= pgprop.Ke;

  return(erate);

}

double get_periglacial_sediment_diffusivity(double Hs,double T0,pgproptype pgprop)
{

  /*uses bilinear interpolation*/

  int i,j;
  double sdiff,ifac,jfac;

  int nHs = pgprop.nHs;
  int nT0 = pgprop.nT0;
  double Kt = pgprop.Kt;

  i=0; j=0;
  while ((pgprop.Hsv[i] < Hs)&&(i < nHs-1)) i += 1;
  while ((pgprop.T0v[j] < T0)&&(j < nT0-1)) j += 1;


  if ((i == 0)||(i == nHs-1)) {
    if ((j == 0)||(j == nT0-1)) sdiff = pgprop.Tr[i][j];
    else { 
      jfac = (T0-pgprop.T0v[j-1])/(pgprop.T0v[j]-pgprop.T0v[j-1]);
      sdiff = (1.0-jfac)*pgprop.Tr[i][j-1] + jfac*pgprop.Tr[i][j];
    }
  }
  else {
    if ((j == 0)||(j == nT0-1)) {
      ifac = (Hs-pgprop.Hsv[i-1])/(pgprop.Hsv[i]-pgprop.Hsv[i-1]);
      sdiff = (1-ifac)*pgprop.Tr[i-1][j]+ifac*pgprop.Tr[i][j];
    }
    else {
      ifac = (Hs-pgprop.Hsv[i-1])/(pgprop.Hsv[i]-pgprop.Hsv[i-1]);
      jfac = (T0-pgprop.T0v[j-1])/(pgprop.T0v[j]-pgprop.T0v[j-1]); 
      sdiff = pgprop.Tr[i-1][j-1]*(1.0-ifac)*(1.0-jfac)+pgprop.Tr[i-1][j]*jfac*(1.0-ifac)+pgprop.Tr[i][j-1]*ifac*(1-jfac)+pgprop.Tr[i][j]*ifac*jfac;
    }
  }
    
  /*scale diffusivity*/
  sdiff *= Kt;
  return(sdiff);

}

void periglacial(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,mproptype mprop,pgproptype pgprop,double dt)
{

  int i,j;
  double sdiff,hdiff,dHs,Ts,kappa,meanerate=0.0,mean_dHs=0.0;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;

  double lrate = mprop.lrate;
  double T0 = mprop.T0;

  double rho_b = pgprop.rho_b;
  double rho_s = pgprop.rho_s;
  double maxice = pgprop.maxice;
  double maxsedi = pgprop.maxsedi;

#pragma omp parallel shared(cells,hp,vp,pgprop,meanerate,mean_dHs) private(i,j,hdiff,sdiff,dHs,Ts,kappa) firstprivate(nx,ny,rho_b,rho_s,dt,lrate,T0,maxsedi,maxice)
  {


    /*initialize*/
    /*#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dHs = 0.0;
	cells[i][j].kappa = 0.0;
      }
      }*/


    /***** loop cells for Erosion *******/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	Ts = T0 - lrate*cells[i][j].topsedi;
	if ((cells[i][j].sedi < maxsedi)&&(cells[i][j].ice < maxice)&&(cells[i][j].fixflag <= 0)) cells[i][j].periglacial_erate = get_periglacial_erosion_rate(cells[i][j].sedi,Ts,cells[i][j].tslope,pgprop);
	else cells[i][j].periglacial_erate = 0.0;
	cells[i][j].periglacial_erosion += cells[i][j].periglacial_erate*dt;
	cells[i][j].bed -= cells[i][j].periglacial_erate*dt;
	cells[i][j].sedi += rho_b/rho_s*cells[i][j].periglacial_erate*dt;
      }
    }
    
    /*loop h-points for horizontal transport*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	Ts = T0 - lrate*(cells[i][j-1].topsedi+cells[i][j].topsedi)/2.0;
	sdiff = get_periglacial_sediment_diffusivity((cells[i][j-1].sedi+cells[i][j].sedi)/2.0,Ts,pgprop);
	dHs = -sdiff*hp[i][j].dbdx*dt/dx;
	if (dHs < 0.0) {
	  if (cells[i][j].sedi <= 0.0) dHs = 0.0;
	  else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	}
	else {
	  if (cells[i][j-1].sedi <= 0.0) dHs = 0.0;
	  else if (dHs > cells[i][j-1].sedi) dHs = cells[i][j-1].sedi;
	}
	hp[i][j].dHs = dHs;
	if (fabs(dHs) > 1.0e-16) hp[i][j].kappa = fabs(dHs)*dx/(fabs(hp[i][j].dbdx)*dt);
	else hp[i][j].kappa = 0.0;
	/*
#pragma omp critical
	{
	  cells[i][j].dHs += dHs;
	  cells[i][j-1].dHs -= dHs;
	  if (fabs(dHs) > 1.0e-16) {
	    kappa = fabs(dHs)*dx/(fabs(hp[i][j].dbdx)*dt);
	    cells[i][j].kappa += sdiff/4.0;
	    cells[i][j-1].kappa += sdiff/4.0;
	  }
	  }*/
      }
    }

    /*loop v-points*/
#pragma omp for schedule(static)
    for (j=0;j<nx;j++) {
      for (i=1;i<ny;i++) {
	Ts = T0 - lrate*(cells[i-1][j].topsedi+cells[i][j].topsedi)/2.0;
	sdiff = get_periglacial_sediment_diffusivity((cells[i-1][j].sedi+cells[i][j].sedi)/2.0,Ts,pgprop);
	dHs = -sdiff*vp[i][j].dbdy*dt/dy;
	if (dHs < 0.0) {
	  if (cells[i][j].sedi <= 0.0) dHs = 0.0; 
	  else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	}
	else {
	  if (cells[i-1][j].sedi <= 0.0) dHs = 0.0;
	  else if (dHs > cells[i-1][j].sedi) dHs = cells[i-1][j].sedi;
	}
	vp[i][j].dHs = dHs;
	if (fabs(dHs) > 1.0e-16) vp[i][j].kappa = fabs(dHs)*dy/(fabs(vp[i][j].dbdy)*dt);
	else vp[i][j].kappa = 0.0;
	/*
#pragma omp critical
	{
	  cells[i][j].dHs += dHs;
	  cells[i-1][j].dHs -= dHs;
	  if (fabs(dHs) > 1.0e-16) {
	    kappa = fabs(dHs)*dy/(fabs(vp[i][j].dbdy)*dt);
	    cells[i][j].kappa += kappa/4.0;
	    cells[i-1][j].kappa += kappa/4.0;
	  }
	  }*/
      }
    }
    
    
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].sedi = cells[i][j].sedi + hp[i][j].dHs - hp[i][j+1].dHs + vp[i][j].dHs - vp[i+1][j].dHs;
	/*cells[i][j].sedi = cells[i][j].sedi + cells[i][j].dHs;*/
	cells[i][j].kappa = 0.25*(hp[i][j].kappa+hp[i][j+1].kappa+vp[i][j].kappa+vp[i+1][j].kappa);
      }
    }
    
  /*compute mean erosion rate and transport change*/
#pragma omp for schedule(static) reduction(+:meanerate) reduction(+:mean_dHs)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanerate += cells[i][j].periglacial_erate;
      mean_dHs += fabs(cells[i][j].kappa);
    }
  }
#pragma omp single
  {
    meanerate /= (double)((*mesh).nc);
    (*mesh).mean_periglacial_erate = meanerate;
    mean_dHs /= (double)((*mesh).nc)*dt*2.0; 
    (*mesh).mean_dHs_periglacial = mean_dHs;
  }


  }

}
