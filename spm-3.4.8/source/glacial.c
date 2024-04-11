

void isosia(celltype **cells,hptype **hp,vptype **vp,cornertype **cp,meshtype *mesh,iproptype iprop)
{

  int i,j;
  long stressit,veloitt;
  double sres,sres_old,sdiff,stot,vres,vres_old,vdiff,vtot,fac;
  double tx,ty,ts,vb,beta,vx,vy,Res,lmaxRes,maxRes = 1.0;
  double dtfac,lmaxdt,maxdt,dum,dhdx,dhdy;
  double maxspeed=0.0,lmaxspeed=0.0;
  double minice = 5.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int periodic = 0;
  int dosliding = (*mesh).dosliding;

  int smode = (*mesh).slidingmode;
  int coldbased = (*mesh).coldbased;
  int gmode = 1; 
  int divcount = 0;

  double gamma = iprop.gamma;
  double gamma0 = iprop.gamma0;
  double ifac = iprop.ifac;
  double sfac = iprop.sfac;
  double vbfac = iprop.vbfac;
  int maxitt_v = iprop.maxitt_v;
  int maxitt_s = iprop.maxitt_s;
  double Cs = iprop.Cs;
  double minbeta = iprop.minbeta;

  double C = iprop.C;
  double L0 = iprop.L0;
  double fc = 1000.0;
  long ndiv = 0;

  double mdh = 0.7;
  double mds = 0.01;


  /*initialize*/
  maxspeed = 0.1;
  lmaxdt = (*mesh).maxdt;
  maxdt = (*mesh).maxdt;



#pragma omp parallel shared(cells,vp,hp,cp,mesh,stressit,sres,sres_old,sdiff,stot,veloitt,vres,vres_old,vdiff,vtot,maxspeed,maxRes,maxdt,ndiv,divcount,ifac) private(i,j,tx,ty,ts,vb,beta,vx,vy,Res,lmaxRes,dtfac,fac,dhdx,dhdy) firstprivate(nx,ny,dx,dy,periodic,lmaxspeed,gamma,gamma0,Cs,sfac,vbfac,maxitt_v,maxitt_s,C,L0,minbeta,dosliding,smode,gmode,minice,lmaxdt,dum,fc,mdh,mds)
  {


#pragma omp single
    {
      veloitt = 0;
      vres = 1.0;
    }
 
    /*itterate velocities*/
      while ((vres > 1.0e-3)&&(veloitt < maxitt_v)) {
      
      /******* Strain rates ********/
      /*shear strain rates*/
#pragma omp for schedule(static)
	for (i=1;i<ny;i++) {
	  for (j=1;j<nx;j++) {
	    cp[i][j].exy = 0.5*((hp[i][j].vx-hp[i-1][j].vx)/dy + (vp[i][j].vy-vp[i][j-1].vy)/dx);
	  } 
	  cp[i][0].exy = 0.5*((hp[i][0].vx-hp[i-1][0].vx)/dy + (vp[i][1].vy-vp[i][0].vy)/dx);
	  cp[i][nx].exy = 0.5*((hp[i][nx].vx-hp[i-1][nx].vx)/dy + (vp[i][nx-1].vy-vp[i][nx-2].vy)/dx);
	}


      /*strain rates in cells*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  if ((cells[i][j].ice > minice)&&(gmode > 0)) {
	    cells[i][j].exx = (hp[i][j+1].vx-hp[i][j].vx)/dx-(cells[i][j].dhdx*((hp[i][j].vx_s-hp[i][j].vx)+(hp[i][j+1].vx_s-hp[i][j+1].vx))/2.0+cells[i][j].dbdx*((hp[i][j].vx-hp[i][j].vx_b)+(hp[i][j+1].vx-hp[i][j+1].vx_b))/2.0)/(cells[i][j].ice+100.0);
	    cells[i][j].eyy = (vp[i+1][j].vy-vp[i][j].vy)/dy-(cells[i][j].dhdy*((vp[i][j].vy_s-vp[i][j].vy)+(vp[i+1][j].vy_s-vp[i+1][j].vy))/2.0+cells[i][j].dbdy*((vp[i][j].vy-vp[i][j].vy_b)+(vp[i+1][j].vy-vp[i+1][j].vy_b))/2.0)/(cells[i][j].ice+100.0);
	    cells[i][j].exy = 0.25*(cp[i][j].exy+cp[i+1][j].exy+cp[i+1][j+1].exy+cp[i][j+1].exy)-(cells[i][j].dhdy*((hp[i][j].vx_s-hp[i][j].vx)+(hp[i][j+1].vx_s-hp[i][j+1].vx))/2.0+cells[i][j].dhdx*((vp[i][j].vy_s-vp[i][j].vy)+(vp[i+1][j].vy_s-vp[i+1][j].vy))/2.0+cells[i][j].dbdy*((hp[i][j].vx-hp[i][j].vx_b)+(hp[i][j+1].vx-hp[i][j+1].vx_b))/2.0+cells[i][j].dbdx*((vp[i][j].vy-vp[i][j].vy_b)+(vp[i+1][j].vy-vp[i+1][j].vy_b))/2.0)/(2.0*cells[i][j].ice+100.0);
	    cells[i][j].ezz = -(cells[i][j].exx+cells[i][j].eyy); 
	  }
	  else {
	    cells[i][j].exx = 0.0;
	    cells[i][j].eyy = 0.0;
	    cells[i][j].exy = 0.0;
	    cells[i][j].ezz = 0.0;
	  }
	}
      }
      
      
#pragma omp single
      {
	stressit = 0;
	sres = 1.0;
      }
      

      /*itterate stress*/
      while ((sres > 1.0e-3)&&(stressit < maxitt_s)) {
	
#pragma omp single
	{
	  sdiff = 0.0;
	  stot = 0.0;
	}
	
	/*compute horizontal stress*/
#pragma omp for schedule(static)
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    if ((cells[i][j].ice > minice)&&(gmode > 0)) {
	      fac = gamma*cells[i][j].te2; 
	      /*if (cells[i][j].ice < 100.0) fac /= cells[i][j].ice/100.0; */
	      if (fac < gamma0) fac = gamma0;
	      cells[i][j].sxx = cells[i][j].exx/fac; if (fabs(cells[i][j].sxx) > fc) cells[i][j].sxx = cells[i][j].sxx/fabs(cells[i][j].sxx)*fc;
	      cells[i][j].syy = cells[i][j].eyy/fac; if (fabs(cells[i][j].syy) > fc) cells[i][j].syy = cells[i][j].syy/fabs(cells[i][j].syy)*fc;
	      cells[i][j].sxy = cells[i][j].exy/fac; if (fabs(cells[i][j].sxy) > fc) cells[i][j].sxy = cells[i][j].sxy/fabs(cells[i][j].sxy)*fc;
	      cells[i][j].szz = -(cells[i][j].sxx+cells[i][j].syy);
	      cells[i][j].td2 = pow(cells[i][j].sxx,2.0)+pow(cells[i][j].syy,2.0)+cells[i][j].sxx*cells[i][j].syy+pow(cells[i][j].sxy,2.0);
	    }
	    else {
	      cells[i][j].sxx = 0.0;
	      cells[i][j].syy = 0.0;
	      cells[i][j].sxy = 0.0;
	      cells[i][j].szz = 0.0;
	      cells[i][j].td2 = 0.0;
	    }
	  } 
	}
	
	/*compute stress gradients - x dir*/
#pragma omp for schedule(static) nowait
	for (i=0;i<ny;i++) {
	  for (j=1;j<nx;j++) {
	    hp[i][j].dsxxdx = (cells[i][j].sxx-cells[i][j-1].sxx)/dx;
	    hp[i][j].dsyydx = (cells[i][j].syy-cells[i][j-1].syy)/dx;
	    hp[i][j].dsxydx = (cells[i][j].sxy-cells[i][j-1].sxy)/dx;
	  }
	  if (hp[i][j].dsxxdx < -mds) hp[i][j].dsxxdx = -mds;  if (hp[i][j].dsxxdx > mds) hp[i][j].dsxxdx = mds; 
	  if (hp[i][j].dsyydx < -mds) hp[i][j].dsyydx = -mds;  if (hp[i][j].dsyydx > mds) hp[i][j].dsyydx = mds; 
	  if (hp[i][j].dsxydx < -mds) hp[i][j].dsxydx = -mds;  if (hp[i][j].dsxydx > mds) hp[i][j].dsxydx = mds; 
	  /*hp[i][0].dsxxdx = hp[i][1].dsxxdx;
	  hp[i][0].dsyydx = hp[i][1].dsyydx;
	  hp[i][0].dsxydx = hp[i][1].dsxydx;
	  hp[i][nx].dsxxdx = hp[i][nx-1].dsxxdx;
	  hp[i][nx].dsyydx = hp[i][nx-1].dsyydx;
	  hp[i][nx].dsxydx = hp[i][nx-1].dsxydx;*/
	  hp[i][0].dsxxdx = 0.0; 
	  hp[i][0].dsyydx = 0.0;
	  hp[i][0].dsxydx = 0.0;
	  hp[i][nx].dsxxdx = 0.0;
	  hp[i][nx].dsyydx = 0.0;
	  hp[i][nx].dsxydx = 0.0;
	}
	
	/*compute stress gradients - y dir*/
#pragma omp for schedule(static)
	for (j=0;j<nx;j++) {
	  for (i=1;i<ny;i++) {
	    vp[i][j].dsxxdy = (cells[i][j].sxx-cells[i-1][j].sxx)/dy;
	    vp[i][j].dsyydy = (cells[i][j].syy-cells[i-1][j].syy)/dy;
	    vp[i][j].dsxydy = (cells[i][j].sxy-cells[i-1][j].sxy)/dy;
	  }
	  if (vp[i][j].dsxxdy < -mds) vp[i][j].dsxxdy = -mds;  if (vp[i][j].dsxxdy > mds) vp[i][j].dsxxdy = mds; 
	  if (vp[i][j].dsyydy < -mds) vp[i][j].dsyydy = -mds;  if (vp[i][j].dsyydy > mds) vp[i][j].dsyydy = mds; 
	  if (vp[i][j].dsxydy < -mds) vp[i][j].dsxydy = -mds;  if (vp[i][j].dsxydy > mds) vp[i][j].dsxydy = mds; 
	  /*vp[0][j].dsxxdy = vp[1][j].dsxxdy;
	  vp[0][j].dsyydy = vp[1][j].dsyydy;
	  vp[0][j].dsxydy = vp[1][j].dsxydy;
	  vp[ny][j].dsxxdy = vp[ny-1][j].dsxxdy;
	  vp[ny][j].dsyydy = vp[ny-1][j].dsyydy;
	  vp[ny][j].dsxydy = vp[ny-1][j].dsxydy;*/
	  vp[0][j].dsxxdy = 0.0;
	  vp[0][j].dsyydy = 0.0;
	  vp[0][j].dsxydy = 0.0;
	  vp[ny][j].dsxxdy = 0.0;
	  vp[ny][j].dsyydy = 0.0;
	  vp[ny][j].dsxydy = 0.0;
	}
	
	/*transfer stress gradients compute coeeficients and update stress*/
#pragma omp for schedule(static)
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    cells[i][j].dsxxdx = 0.5*(hp[i][j].dsxxdx+hp[i][j+1].dsxxdx);
	    cells[i][j].dsyydx = 0.5*(hp[i][j].dsyydx+hp[i][j+1].dsyydx);
	    cells[i][j].dsxydx = 0.5*(hp[i][j].dsxydx+hp[i][j+1].dsxydx);
	  }
	}
#pragma omp for schedule(static)
	for (j=0;j<nx;j++) {
	  for (i=0;i<ny;i++) {
	    cells[i][j].dsxxdy = 0.5*(vp[i][j].dsxxdy+vp[i+1][j].dsxxdy);
	    cells[i][j].dsyydy = 0.5*(vp[i][j].dsyydy+vp[i+1][j].dsyydy);
	    cells[i][j].dsxydy = 0.5*(vp[i][j].dsxydy+vp[i+1][j].dsxydy);
	  }
	}


#pragma omp for schedule(static) reduction(+:sdiff) reduction(+:stot)
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {	 
	    cells[i][j].cx[0] = (2.0*cells[i][j].sxx+cells[i][j].syy)*cells[i][j].dhdx+cells[i][j].sxy*cells[i][j].dhdy;
	    cells[i][j].cx[1] = -(1.0+cells[i][j].alpha2)*cells[i][j].dhdx+2.0*cells[i][j].dsxxdx+cells[i][j].dsyydx+cells[i][j].dsxydy;
	    cells[i][j].cy[0] = (cells[i][j].sxx+2.0*cells[i][j].syy)*cells[i][j].dhdy+cells[i][j].sxy*cells[i][j].dhdx;
	    cells[i][j].cy[1] = -(1.0+cells[i][j].alpha2)*cells[i][j].dhdy+cells[i][j].dsxxdy+2.0*cells[i][j].dsyydy+cells[i][j].dsxydx;
	    cells[i][j].kc[0] = cells[i][j].cx[0]*cells[i][j].cx[0]+cells[i][j].cy[0]*cells[i][j].cy[0]+cells[i][j].td2;
	    cells[i][j].kc[1] = 2.0*(cells[i][j].cx[0]*cells[i][j].cx[1]+cells[i][j].cy[0]*cells[i][j].cy[1]);
	    cells[i][j].kc[2] = cells[i][j].cx[1]*cells[i][j].cx[1]+cells[i][j].cy[1]*cells[i][j].cy[1];
	    cells[i][j].te2r = cells[i][j].kc[0]+cells[i][j].kc[1]*cells[i][j].ice/2.0+cells[i][j].kc[2]*pow(cells[i][j].ice,2.0)/3.0;
	    cells[i][j].te2 += sfac*(cells[i][j].te2r - cells[i][j].te2);
	    sdiff += pow(cells[i][j].te2r - cells[i][j].te2,2.0);
	    stot += pow(cells[i][j].te2r,2.0);
	  }
	}
	
	
	  
#pragma omp single 
	{
	  sres_old = sres;
	  sres = sdiff/(stot+1.0e-16);
	  stressit += 1;
	}

      }/*while*/      
      
#pragma omp single 
      {
	(*mesh).stressitt += stressit;
	(*mesh).stressitt_count += 1;
      }


      /************* Velocity update ***************/
	
#pragma omp single
      {
	vdiff = 0.0;
	vtot = 0.0;
	maxRes = 0.0;
      }
      
      lmaxRes = 0.0;

      /*x velocity*/
#pragma omp for schedule(static) reduction(+:vdiff) reduction(+:vtot)
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) {
	  if (hp[i][j].ice > minice) {
	    dhdx = hp[i][j].dhdx; if (dhdx < -mdh) dhdx = -mdh; if (dhdx > mdh) dhdx = mdh; 
	    dhdy = hp[i][j].dhdy; if (dhdy < -mdh) dhdy = -mdh; if (dhdy > mdh) dhdy = mdh; 
	    hp[i][j].cx[0] = (cells[i][j-1].sxx+cells[i][j].sxx+.5*(cells[i][j-1].syy+cells[i][j].syy))*dhdx+.5*(cells[i][j-1].sxy+cells[i][j].sxy)*dhdy;
	    hp[i][j].cx[1] = -(1.0+hp[i][j].alpha2)*dhdx+2.0*hp[i][j].dsxxdx+hp[i][j].dsyydx+.5*(cells[i][j-1].dsxydy+cells[i][j].dsxydy); 
	    hp[i][j].cy[0] = (.5*(cells[i][j-1].sxx+cells[i][j].sxx)+cells[i][j-1].syy+cells[i][j].syy)*dhdy+.5*(cells[i][j-1].sxy+cells[i][j].sxy)*dhdx;
	    hp[i][j].cy[1] = -(1.0+hp[i][j].alpha2)*dhdy+.5*(cells[i][j-1].dsxxdy+cells[i][j].dsxxdy)+cells[i][j-1].dsyydy+cells[i][j].dsyydy+hp[i][j].dsxydx;
	    hp[i][j].kc[0] = hp[i][j].cx[0]*hp[i][j].cx[0]+hp[i][j].cy[0]*hp[i][j].cy[0]+.5*(cells[i][j-1].td2+cells[i][j].td2);
	    hp[i][j].kc[1] = 2.0*(hp[i][j].cx[0]*hp[i][j].cx[1]+hp[i][j].cy[0]*hp[i][j].cy[1]);
	    hp[i][j].kc[2] = hp[i][j].cx[1]*hp[i][j].cx[1]+hp[i][j].cy[1]*hp[i][j].cy[1];
	    hp[i][j].wx[0] = hp[i][j].cx[0]*hp[i][j].kc[0];/*-((hp[i][j+1].vz_b-hp[i][j-1].vz_b)/(2.0*dx)-hp[i][j].dbdx*(cells[i][j-1].ezz+cells[i][j].ezz)/2.0)/(2.0*gamma);*/
	    hp[i][j].wx[1] = hp[i][j].cx[0]*hp[i][j].kc[1]+hp[i][j].cx[1]*hp[i][j].kc[0];/*-1.0/(4.0*gamma)*(cells[i][j].ezz-cells[i][j-1].ezz)/dx;*/
	    hp[i][j].wx[2] = hp[i][j].cx[0]*hp[i][j].kc[2]+hp[i][j].cx[1]*hp[i][j].kc[1];
	    hp[i][j].wx[3] = hp[i][j].cx[1]*hp[i][j].kc[2]; 
	    hp[i][j].vrx = 2.0*gamma*(hp[i][j].wx[0]*hp[i][j].ice/2.0+hp[i][j].wx[1]*pow(hp[i][j].ice,2.0)/3.0+hp[i][j].wx[2]*pow(hp[i][j].ice,3.0)/4.0+hp[i][j].wx[3]*pow(hp[i][j].ice,4.0)/5.0);
	    hp[i][j].dvx_d = ifac*(hp[i][j].vrx - hp[i][j].vx_d);
	    hp[i][j].vx_d += hp[i][j].dvx_d;
	    Res = fabs(hp[i][j].vrx-hp[i][j].vx_d)/(fabs(hp[i][j].vx_d)+1.0);
	    if (Res > lmaxRes) lmaxRes = Res;
	    hp[i][j].vresx = hp[i][j].vrx - hp[i][j].vx_d;
	    if ((cells[i][j-1].margin < 0)&&(cells[i][j].margin < 0)) {
	      vdiff += pow(hp[i][j].vrx-hp[i][j].vx_d,2.0);
	      vtot += pow(hp[i][j].vrx,2.0); 
	    }
	  }
	  else {
	    hp[i][j].vx_d = 0.0;
	    hp[i][j].dvx_d = 0.0;
	    hp[i][j].vresx = 0.0;;
	  }
	}

	hp[i][0].vx_d = 0.0;
	hp[i][nx].vx_d = 0.0;

	/*hp[i][0].vx_d = 2.0*hp[i][1].vx_d-hp[i][2].vx_d;
	hp[i][nx].vx_d = 2.0*hp[i][nx-1].vx_d-hp[i][nx-2].vx_d;
	if ((hp[i][0].vx_d*hp[i][1].vx_d) < 0.0) hp[i][0].vx_d = 0.0;
	if ((hp[i][nx].vx_d*hp[i][nx-1].vx_d) < 0.0) hp[i][nx].vx_d = 0.0;*/
      }



      
      /*y velocity*/
#pragma omp for schedule(static) reduction(+:vdiff) reduction(+:vtot)
      for (j=0;j<nx;j++) {
	for (i=1;i<ny;i++) {
	  if (vp[i][j].ice > minice) {
	    dhdx = vp[i][j].dhdx; if (dhdx < -mdh) dhdx = -mdh; if (dhdx > mdh) dhdx = mdh; 
	    dhdy = vp[i][j].dhdy; if (dhdy < -mdh) dhdy = -mdh; if (dhdy > mdh) dhdy = mdh; 
	    vp[i][j].cx[0] = (cells[i-1][j].sxx+cells[i][j].sxx+.5*(cells[i-1][j].syy+cells[i][j].syy))*dhdx+.5*(cells[i-1][j].sxy+cells[i][j].sxy)*dhdy;
	    vp[i][j].cx[1] = -(1.0+vp[i][j].alpha2)*dhdx+(cells[i-1][j].dsxxdx+cells[i][j].dsxxdx)+.5*(cells[i-1][j].dsyydx+cells[i][j].dsyydx)+vp[i][j].dsxydy;
	    vp[i][j].cy[0] = (.5*(cells[i-1][j].sxx+cells[i][j].sxx)+cells[i-1][j].syy+cells[i][j].syy)*dhdy+.5*(cells[i-1][j].sxy+cells[i][j].sxy)*dhdx;
	    vp[i][j].cy[1] = -(1.0+vp[i][j].alpha2)*dhdy+vp[i][j].dsxxdy+2.0*vp[i][j].dsyydy+.5*(cells[i-1][j].dsxydx+cells[i][j].dsxydx);
	    vp[i][j].kc[0] = vp[i][j].cx[0]*vp[i][j].cx[0]+vp[i][j].cy[0]*vp[i][j].cy[0]+.5*(cells[i-1][j].td2+cells[i][j].td2);
	    vp[i][j].kc[1] = 2.0*(vp[i][j].cx[0]*vp[i][j].cx[1]+vp[i][j].cy[0]*vp[i][j].cy[1]);
	    vp[i][j].kc[2] = vp[i][j].cx[1]*vp[i][j].cx[1]+vp[i][j].cy[1]*vp[i][j].cy[1];
	    vp[i][j].wy[0] = vp[i][j].cy[0]*vp[i][j].kc[0];/*-((vp[i+1][j].vz_b-vp[i-1][j].vz_b)/(2.0*dy)-vp[i][j].dbdy*(cells[i-1][j].ezz+cells[i][j].ezz)/2.0)/(2.0*gamma);*/
	    vp[i][j].wy[1] = vp[i][j].cy[0]*vp[i][j].kc[1]+vp[i][j].cy[1]*vp[i][j].kc[0];/*-1.0/(4.0*gamma)*(cells[i][j].ezz-cells[i-1][j].ezz)/dy;*/
	    vp[i][j].wy[2] = vp[i][j].cy[0]*vp[i][j].kc[2]+vp[i][j].cy[1]*vp[i][j].kc[1];
	    vp[i][j].wy[3] = vp[i][j].cy[1]*vp[i][j].kc[2];
	    vp[i][j].vry = 2.0*gamma*(vp[i][j].wy[0]*vp[i][j].ice/2.0+vp[i][j].wy[1]*pow(vp[i][j].ice,2.0)/3.0+vp[i][j].wy[2]*pow(vp[i][j].ice,3.0)/4.0+vp[i][j].wy[3]*pow(vp[i][j].ice,4.0)/5.0);
	    vp[i][j].dvy_d = ifac*(vp[i][j].vry - vp[i][j].vy_d);
	    vp[i][j].vy_d += vp[i][j].dvy_d;
	    Res = fabs(vp[i][j].vry-vp[i][j].vy_d)/(fabs(vp[i][j].vy_d)+1.0);
	    if (Res > lmaxRes) lmaxRes = Res;
	    vp[i][j].vresy = vp[i][j].vry - vp[i][j].vy_d; 
	    if ((cells[i-1][j].margin < 0)&&(cells[i][j].margin < 0)) {
	      vdiff += pow(vp[i][j].vry-vp[i][j].vy_d,2.0);
	      vtot += pow(vp[i][j].vry,2.0); 
	    }
	  }
	  else {
	    vp[i][j].vy_d = 0.0;
	    vp[i][j].dvy_d = 0.0;
	    vp[i][j].vresy = 0.0;
	  }
	}
	vp[0][j].vy_d = 0.0;/*vp[1][j].vy_d;*/
	vp[ny][j].vy_d = 0.0;/*vp[ny-1][j].vy_d;*/
      } 


#pragma omp critical 
      { 
	if (lmaxRes > maxRes) maxRes = lmaxRes; 
      }

      /******** basal sliding **********/
      if (dosliding > 0) {

#pragma omp for schedule(static) 
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) { 
	    if (cells[i][j].margin < 0.0) {
	      cells[i][j].lb = sqrt(1.0+pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0));
	      /*cells[i][j].pb = (1.0+pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0))*cells[i][j].ice+cells[i][j].szz;*/
	      cells[i][j].pb = cells[i][j].ice+cells[i][j].szz;
	    cells[i][j].sxz = cells[i][j].cx[0]+cells[i][j].cx[1]*cells[i][j].ice;
	    cells[i][j].syz = cells[i][j].cy[0]+cells[i][j].cy[1]*cells[i][j].ice;
	    cells[i][j].tn = cells[i][j].pb+(2.0*cells[i][j].dbdx*cells[i][j].sxz+2.0*cells[i][j].dbdy*cells[i][j].syz-pow(cells[i][j].dbdx,2.0)*cells[i][j].sxx-pow(cells[i][j].dbdy,2.0)*cells[i][j].syy-cells[i][j].szz-2.0*cells[i][j].dbdx*cells[i][j].dbdy*cells[i][j].sxy)/pow(cells[i][j].lb,2.0);
	    /*cells[i][j].tn = cells[i][j].ice;*/
	    if (cells[i][j].tn < 0.0) cells[i][j].tn = 0.0;
	    /*cells[i][j].te = cells[i][j].tn-cells[i][j].Pw; if (cells[i][j].te < 0.01*cells[i][j].ice) cells[i][j].te = 0.01*cells[i][j].ice;*/
	    cells[i][j].tbx = (cells[i][j].dbdx*(cells[i][j].pb-cells[i][j].sxx)-cells[i][j].dbdy*cells[i][j].sxy+cells[i][j].sxz-cells[i][j].dbdx*cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].tby = (-cells[i][j].dbdx*cells[i][j].sxy+cells[i][j].dbdy*(cells[i][j].pb-cells[i][j].syy)+cells[i][j].syz-cells[i][j].dbdy*cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].tbz = (-cells[i][j].dbdx*cells[i][j].sxz-cells[i][j].dbdy*cells[i][j].syz-cells[i][j].pb+cells[i][j].szz+cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].ts = sqrt(pow(cells[i][j].tbx,2.0)+pow(cells[i][j].tby,2.0)+pow(cells[i][j].tbz,2.0));
	    if (smode == 0) {
	      if (cells[i][j].ice > 10.0) cells[i][j].vb = cells[i][j].Cs*pow(cells[i][j].ts,2.0);
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (cells[i][j].Cs*cells[i][j].ice/10.0)*pow(cells[i][j].ts,2.0);
	      else cells[i][j].vb = 0.0;
	    }
	    else if (smode == 1) {
	      if (cells[i][j].ice > 10.0) cells[i][j].vb = cells[i][j].Cs*pow(cells[i][j].ts,2.0)/(cells[i][j].te+1.0);
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (cells[i][j].Cs*cells[i][j].ice/10.0)*pow(cells[i][j].ts,2.0)/(cells[i][j].te+1.0);
	      else cells[i][j].vb = 0.0;
	    }
	    else if (smode == 2) {
	      beta = pow(C,3.0)-pow(cells[i][j].ts/(cells[i][j].te+1.0),3.0);
	      if (beta < minbeta) beta = minbeta;
	      else if (cells[i][j].ice > 10.0) cells[i][j].vb = L0*pow(cells[i][j].ts,3.0)/beta;
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (L0*cells[i][j].ice/10.0)*pow(cells[i][j].ts,3.0)/beta;
	      else cells[i][j].vb = 0.0;	    
	    }
	    else if (smode == 3) {
	      beta = 1.0-cells[i][j].SLf;
	      if (beta < 0.1) beta = 0.1;
	      if (cells[i][j].ice > 100.0) cells[i][j].vb = cells[i][j].Cs*pow(cells[i][j].ts,3.0)/beta;
	      else if (cells[i][j].ice > 10.0) cells[i][j].vb = (cells[i][j].Cs*cells[i][j].ice/100.0)*pow(cells[i][j].ts,3.0)/beta;
	      else cells[i][j].vb = 0.0; 
	    }	      
	    else {
	      dum = cells[i][j].ts+cells[i][j].Pw*cells[i][j].slidingslope;
	      if (dum < 0.1*cells[i][j].ts) dum = 0.1*cells[i][j].ts;
	      beta = 1.0-cells[i][j].SLf;
	      if (beta < 1.0e-2) beta = 1.0e-2;
	      if (cells[i][j].ice > 10.0) cells[i][j].vb = cells[i][j].Cs*pow(dum,3.0)/beta;
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (cells[i][j].Cs*cells[i][j].ice/10.0)*pow(dum,3.0)/beta;
	      else cells[i][j].vb = 0.0;
	    }
	    if (coldbased == 1) cells[i][j].vb *= cells[i][j].sfac; 
	     } 
	       else cells[i][j].vb = 0.0;
	  }
	}

	/*transfer basal sliding*/
#pragma omp for schedule(static) nowait
	for (i=0;i<ny;i++) {
	  for (j=1;j<nx;j++) {
	    tx = .5*(cells[i][j-1].tbx+cells[i][j].tbx);
	    ts = .5*(cells[i][j-1].ts+cells[i][j].ts);
	    vb = .5*(cells[i][j-1].vb+cells[i][j].vb);
	    if (hp[i][j].ice > 1.0) hp[i][j].vx_b = (1.0-vbfac)*hp[i][j].vx_b + vbfac*vb*tx/(ts+1.0e-6);
	    else hp[i][j].vx_b = 0.0;
	  }
	  hp[i][0].vx_b = 2.0*hp[i][1].vx_b-hp[i][2].vx_b;
	  hp[i][nx].vx_b = 2.0*hp[i][nx-1].vx_b-hp[i][nx-2].vx_b;
	  if ((hp[i][0].vx_b*hp[i][1].vx_b) < 0.0) hp[i][0].vx_b = 0.0;
	  if ((hp[i][nx].vx_b*hp[i][nx-1].vx_b) < 0.0) hp[i][nx].vx_b = 0.0;
	}
#pragma omp for schedule(static)
      for (j=0;j<nx;j++) {
 	for (i=1;i<ny;i++) {
	    ty = .5*(cells[i-1][j].tby+cells[i][j].tby);
	    ts = .5*(cells[i-1][j].ts+cells[i][j].ts);
	    vb = .5*(cells[i-1][j].vb+cells[i][j].vb);
	    if (vp[i][j].ice > 1.0) vp[i][j].vy_b = (1.0-vbfac)*vp[i][j].vy_b + vbfac*vb*ty/(ts+1.0e-6);
	    else vp[i][j].vy_b = 0.0;
	}
	vp[0][j].vy_b = 0.0;/*vp[1][j].vy_b;*/
	vp[ny][j].vy_b = 0.0;/*vp[ny-1][j].vy_b;*/
      } 
      
      }/*dosliding*/

      /*combine velocities*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx+1;j++) {
	  hp[i][j].vx = hp[i][j].vx_d + hp[i][j].vx_b;
	}
      }
#pragma omp for schedule(static)
      for (i=0;i<ny+1;i++) {
	for (j=0;j<nx;j++) {
	  vp[i][j].vy = vp[i][j].vy_d + vp[i][j].vy_b;
	}
      }
      
      /*approximate surface velocities*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx+1;j++) {
	  hp[i][j].vx_s = 1.2*hp[i][j].vx;
	}
      }
#pragma omp for schedule(static)
      for (i=0;i<ny+1;i++) {
	for (j=0;j<nx;j++) {
	  vp[i][j].vy_s = 1.2*vp[i][j].vy; 
	}
      }
      
#pragma omp single
      {
	/*if (vtot < (double)(nx*ny)) vtot = (double)(nx*ny);*/
	vres_old = vres;
	vres = vdiff/(vtot+1e-16);
	veloitt += 1;
	if (vres > vres_old) {
	  ndiv += 1;
	  /*veloitt = maxitt_v + 1;*/
	}
	
      }

    }/*while*/
    
    
    /*compute surface velocities*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx+1;j++) {
	hp[i][j].vx_s = 2.0*gamma*(hp[i][j].wx[0]*hp[i][j].ice+hp[i][j].wx[1]*pow(hp[i][j].ice,2.0)/2.0+hp[i][j].wx[2]*pow(hp[i][j].ice,3.0)/3.0+hp[i][j].wx[3]*pow(hp[i][j].ice,4.0)/4.0);
      }
    }
#pragma omp for schedule(static)
    for (i=0;i<ny+1;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].vy_s = 2.0*gamma*(vp[i][j].wy[0]*vp[i][j].ice+vp[i][j].wy[1]*pow(vp[i][j].ice,2.0)/2.0+vp[i][j].wy[2]*pow(vp[i][j].ice,3.0)/3.0+vp[i][j].wy[3]*pow(vp[i][j].ice,4.0)/4.0); 
      }
    }


  /*transfer velocity info to cells*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].vx_d = .5*(hp[i][j].vx_d+hp[i][j+1].vx_d);
	cells[i][j].vy_d = .5*(vp[i][j].vy_d+vp[i+1][j].vy_d);
	cells[i][j].deformation = sqrt(pow(cells[i][j].vx_d,2.0)+pow(cells[i][j].vy_d,2.0));
	cells[i][j].vx_b = .5*(hp[i][j].vx_b+hp[i][j+1].vx_b);
	cells[i][j].vy_b = .5*(vp[i][j].vy_b+vp[i+1][j].vy_b);
	cells[i][j].sliding = sqrt(pow(cells[i][j].vx_b,2.0)+pow(cells[i][j].vy_b,2.0));
	cells[i][j].vbres = cells[i][j].vb - cells[i][j].sliding;
	cells[i][j].ezz_s = -((hp[i][j+1].vx-hp[i][j].vx)/dx+(vp[i+1][j].vy-vp[i][j].vy)/dy);
      }
    }

    /*find max speed for time step scaling*/
    lmaxspeed = 0.0;
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<=nx;j++) {
	if (fabs(hp[i][j].vx) > lmaxspeed) lmaxspeed = fabs(hp[i][j].vx);
	dtfac = (*mesh).ct*dx*dx*fabs(hp[i][j].dhdx/(hp[i][j].vx*hp[i][j].ice));
	if (lmaxdt > dtfac) lmaxdt = dtfac;
      }
    }
#pragma omp for schedule(static)
    for (i=0;i<=ny;i++) {
      for (j=0;j<nx;j++) {
	if (fabs(vp[i][j].vy) > lmaxspeed) lmaxspeed = fabs(vp[i][j].vy);
	dtfac = (*mesh).ct*dy*dy*fabs(vp[i][j].dhdy/(vp[i][j].vy*vp[i][j].ice));
	if (lmaxdt > dtfac) lmaxdt = dtfac;
      }
    }
#pragma omp critical
    {
      if (lmaxspeed > maxspeed) maxspeed = lmaxspeed;
      if (maxdt > lmaxdt) maxdt = lmaxdt;
    }
    
  }/*pragma*/

  /*printf(" veloitt = %ld ",veloitt);*/
  (*mesh).veloitt += veloitt;
  (*mesh).veloitt_count += 1;
  (*mesh).maxspeed = maxspeed;  
  (*mesh).maxdt_ice = maxdt;

}/*isosia*/

void depositsnow(celltype **cells,int i,int j,double sc,double cc)
{

  int k,nne,maxk;
  double maxslope,curv,maxh,gift;
  double neslope[8];
  double nefac[8];
 
  /*number of neighbours*/
  nne = cells[i][j].nne;

  /*ice surface curvature*/
  curv = cells[i][j].curv;

  /*loop neighbors and find steepest downslope neighbour*/
  neslope[0] = (cells[i][j].topsnow-cells[cells[i][j].ne_i[0]][cells[i][j].ne_j[0]].topsnow)/cells[i][j].ndist[0];
  maxk = 0; maxslope = neslope[0];
  for (k=1;k<nne;k++) {
    neslope[k] = (cells[i][j].topsnow-cells[cells[i][j].ne_i[k]][cells[i][j].ne_j[k]].topsnow)/cells[i][j].ndist[k];
      if (neslope[k] > maxslope) {
	maxslope = neslope[k];
	maxk = k;
      }
    }

    /*if conditions for for avalancing*/
    if ((maxslope > sc)||(curv < cc)) {

      cells[i][j].avasite = 1;

      /*if reciever is not already a doner - avoiding endless loops*/ 
      if (cells[cells[i][j].ne_i[maxk]][cells[i][j].ne_j[maxk]].avasite == 0) {
      
	/*max alowed elevation*/
	if (curv < cc)
	  maxh = cells[i][j].topice;
	else
	  maxh = cells[cells[i][j].ne_i[maxk]][cells[i][j].ne_j[maxk]].topsnow + sc*cells[i][j].ndist[maxk];

	/*snow package to send*/
	gift = cells[i][j].topsnow - maxh;
	if (gift > cells[i][j].snow) gift = cells[i][j].snow;
	if (gift < 0.0) gift = 0.0;

	cells[i][j].snow -= gift;
	cells[i][j].topsnow -= gift;

	cells[cells[i][j].ne_i[maxk]][cells[i][j].ne_j[maxk]].snow += gift;
	cells[cells[i][j].ne_i[maxk]][cells[i][j].ne_j[maxk]].topsnow += gift;

	/*recursive call*/
	depositsnow(cells,cells[i][j].ne_i[maxk],cells[i][j].ne_j[maxk],sc,cc);

      }
      else {
	cells[i][j].snow = 0.0;
	cells[i][j].topsnow = cells[i][j].topice;
      }

    }

}



void avalance(celltype **cells,meshtype mesh,mproptype mprop,sortarraytype *sortarray)
{


  int i,j,k,indx;
  
  int nx = mesh.nx;
  int ny = mesh.ny;
  int nn = nx*ny;
  double dts = 0.1; /*snow accumulation time*/
  double pos = 0.15; /*snow pack density ratio*/ 

  double sc = mprop.avaslope;
  double cc = mprop.avacurv;

  
  /*sort cells in order of decreasing elevation*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].snow = dts*cells[i][j].accrate/pos;
      cells[i][j].topsnow = cells[i][j].topice+cells[i][j].snow;
      cells[i][j].avasite = 0;
      indx = i*nx + j;
      sortarray[indx].ii = i;
      sortarray[indx].jj = j;
      sortarray[indx].value = cells[i][j].topsnow;
    }
  }
  qsort(sortarray,nn,sizeof(sortarraytype),compare_values);

  /*loop cells in decending order of elevation*/ 
  for (indx=0;indx<nn;indx++) {
    i = sortarray[indx].ii;
    j = sortarray[indx].jj;

    /*check for deposit conditions*/
    if (cells[i][j].snow > 1e-3) depositsnow(cells,i,j,sc,cc);

  }/*indx*/

  /*update accumulation rates*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].accrate = pos*cells[i][j].snow/dts;
    }
  }

}


void get_mrate_h(celltype **cells,meshtype mesh,mproptype mprop)
{

  /*interpolate mass balance function based on elevation*/
  int i,j,k;
  double fac,mrate,h;

  int nx = mesh.nx;
  int ny = mesh.ny;

  double maxabla = 10.0;
  
  k = 0;

  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      h = cells[i][j].topice;

      if (h > mprop.Mrate_h[0][k]) {
	while ((h > mprop.Mrate_h[0][k])&&(k < mprop.nMrate_h-1)) k += 1;
      }
      else if (k > 0) {
	while ((h < mprop.Mrate_h[0][k-1])&&(k > 0)) k -= 1;
      }
	
      if ((k == 0)||(k == mprop.nMrate_h-1)) mrate = mprop.Mrate_h[1][k];
      else {
	fac = (h - mprop.Mrate_h[0][k-1])/(mprop.Mrate_h[0][k]-mprop.Mrate_h[0][k-1]);
	mrate = (1.0-fac)*mprop.Mrate_h[1][k-1] + fac*mprop.Mrate_h[1][k];
      }

      if (cells[i][j].include > 0) {
	if (mrate > 0.0) {
	  cells[i][j].accrate = mrate;
	  cells[i][j].smelt = 0.0;
	}
	else {
	  cells[i][j].accrate = 0.0;
	  cells[i][j].smelt = -mrate;
	}
      }
      else {
	cells[i][j].accrate = 0.0;
	cells[i][j].smelt = maxabla;
      }

    }
  }

}

void get_mrate_T(celltype **cells,meshtype mesh,mproptype mprop)
{

  /*interpolate mass balance function based on temperature*/
  int i,j,k;
  double fac,smelt,accrate,T;

  int nx = mesh.nx;
  int ny = mesh.ny;

  double maxabla = 10.0;
  
  k = 0;

  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      T = cells[i][j].Tair;

      if (T > mprop.Mrate_T[0][k]) {
	while ((T > mprop.Mrate_T[0][k])&&(k < mprop.nMrate_T-1)) k += 1;
      }
      else if (k > 0) {
	while ((T < mprop.Mrate_T[0][k-1])&&(k > 0)) k -= 1;
      }
	
      if ((k == 0)||(k == mprop.nMrate_T-1)) {
	accrate = mprop.Mrate_T[1][k];
	smelt = mprop.Mrate_T[2][k];
      }
      else {
	fac = (T - mprop.Mrate_T[0][k-1])/(mprop.Mrate_T[0][k]-mprop.Mrate_T[0][k-1]);
	accrate = (1.0-fac)*mprop.Mrate_T[1][k-1] + fac*mprop.Mrate_T[1][k];
	smelt = (1.0-fac)*mprop.Mrate_T[2][k-1] + fac*mprop.Mrate_T[2][k];
      }

      if (cells[i][j].include > 0) {
	  cells[i][j].accrate = cells[i][j].precip*accrate;
	  cells[i][j].smelt = smelt;
      }
      else {
	cells[i][j].accrate = 0.0;
	cells[i][j].smelt = maxabla;
      }

    }
  }

}

void accumulation_and_melt(celltype **cells,meshtype mesh,mproptype mprop)
{

  int i,j;
  
  int nx = mesh.nx;
  int ny = mesh.ny;
  double mrate;

  int mtype = mprop.mtype;

  /*when mass balance is function of elevation*/
  if (mtype == 1) get_mrate_h(cells,mesh,mprop);

  /*when mass balance is function of temperature*/
  else if (mtype == 2) get_mrate_T(cells,mesh,mprop);

  /*when mass balance is fixed for each cell*/
  else {

    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	mrate = cells[i][j].mrate;
	if (mrate > 0.0) {
	  cells[i][j].accrate = mrate;
	  cells[i][j].smelt = 0.0;
	}
	else {
	  cells[i][j].accrate = 0.0;
	  cells[i][j].smelt = -mrate;
	}

      }/*j*/
    }/*i*/

  }

}

void mass_balance(celltype **cells,meshtype *mesh,iproptype iprop,mproptype mprop,hwproptype hwprop,double time,double dt)
{

  int i,j;
  double melt,accumulation,smelt,bmelt,imelt;
  double meanice_b,dice,hs = 0.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int doglacialsedi = (*mesh).doglacialsedi;
  int dodebrisablation = (*mesh).dodebrisablation;
  double g = (*mesh).gravity;
  double L = iprop.latentheat; 
  double Ldebris = mprop.Ldebris;
  double maxacc = mprop.maxacc;

#pragma omp parallel shared(cells,mesh,mprop,hwprop,dice,meanice_b) private(i,j,melt,smelt,bmelt,imelt,accumulation,hs) firstprivate(nx,ny,time,dt,g,L,doglacialsedi,dodebrisablation,Ldebris,maxacc)
  {

    /*compute acuumulation and ablation*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	
	/*ice accumulation*/
	if (cells[i][j].accrate > maxacc) cells[i][j].accrate = maxacc;
	accumulation = dt*cells[i][j].accrate;

	/*ablation rate*/
	smelt = cells[i][j].smelt;

	/*if (cells[i][j].x < 50.0e3) {
	  smelt = 10.0; 
	}
	else if (cells[i][j].bed < -0.9*cells[i][j].ice) {
	  smelt = 10.0; 
	}
	else {
	  smelt = cells[i][j].smelt;
	}*/

	/*basal melting*/
	bmelt = -cells[i][j].Mb; /*Mb is negative in thermal.c*/

	/*internal melting*/
	imelt = 9.82*sqrt(cells[i][j].te2)*cells[i][j].deformation/L;

	/*reduce surface melt due to debris cover*/
	hs = cells[i][j].ssedi; 
	if (hs > 2.0) hs = 2.0;
	if (hs < 0.0) hs = 0.0;

	/*if (dodebrisablation) smelt *= exp(-cells[i][j].Vs[0]/Ldebris);*/
	if (dodebrisablation) smelt *= Ldebris/(Ldebris+hs);
	
	/*scale melt*/
	smelt *= dt;
	bmelt *= dt; 
	imelt *= dt;

	/*total melt*/
	melt = smelt + bmelt + imelt;

	/*limit melting*/
	if ((melt > (cells[i][j].ice+accumulation))&&(melt > 0.0)) {
	  bmelt *= (cells[i][j].ice+accumulation)/melt;
	  smelt *= (cells[i][j].ice+accumulation)/melt;
	  imelt *= (cells[i][j].ice+accumulation)/melt;
	  melt = cells[i][j].ice + accumulation;
	}

	/*store melt rate*/
	cells[i][j].meltrate = melt/dt;
	cells[i][j].Ms = (accumulation - melt)/dt;

	/*do mass balance*/
	cells[i][j].ice += accumulation - melt;

	/*pass melts to hydrology*/
	if ((*mesh).doglacialhydrology > 0) cells[i][j].water = (hwprop.a2w*smelt + bmelt + imelt)/dt; /*m yr-1*/

      }/*j*/
    }/*i*/


    /*compute total ice volume change*/
#pragma omp for schedule(static) reduction(+:dice)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      dice += cells[i][j].Ms;
    }
  }
#pragma omp single
  {
    meanice_b = dice*dt/((double)(nx*ny));
    (*mesh).meanice_b += meanice_b;
    }

  }/*pragma*/

}

void get_change(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,double dt)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double dH,hdiff;
  double meanice = 0.0;

#pragma omp parallel shared(cells,vp,hp,mesh) private(i,j,dH,hdiff) firstprivate(nx,ny,dx,dy,dt)
  {
    
    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	/*cells[i][j].dH = 0.0;*/
	cells[i][j].oldice = cells[i][j].ice;
      }
    }


  /*loop h-points*/
#pragma omp for schedule(static)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx+1;j++) {
      dH = hp[i][j].vx*hp[i][j].ice*dt/dx;
      if ((dH < 0.0)&&(j<nx)) {
	if (cells[i][j].ice <= 0.0) dH = 0.0;
	else if (dH < -0.25*cells[i][j].ice) dH = -0.25*cells[i][j].ice;
      }
      else if ((dH > 0.0)&&(j > 0)) {
	if (cells[i][j-1].ice <= 0.0) dH = 0.0;
	else if (dH > 0.25*cells[i][j-1].ice) dH = 0.25*cells[i][j-1].ice;
      }
      hp[i][j].dH = dH;
    }
  }

  /*loop v-points*/
#pragma omp for schedule(static)
  for (i=1;i<ny;i++) {
    for (j=0;j<nx;j++) {
      dH = vp[i][j].vy*vp[i][j].ice*dt/dy;
      /*hdiff = .5*(cells[i-1][j].topice-cells[i][j].topice);
      if ((dH > 0.0)&&(dH > hdiff)) { dH = hdiff; if (dH < 0.0) dH = 0.0;}
      if ((dH < 0.0)&&(dH < hdiff)) { dH = hdiff; if (dH > 0.0) dH = 0.0;}
      if ((dH > 0.0)&&(cells[i-1][j].ice < 0.0)) dH = 0.0;
      if ((dH < 0.0)&&(cells[i][j].ice < 0.0)) dH = 0.0;*/
      if (dH < 0.0) {
	if (cells[i][j].ice <= 0.0) dH = 0.0; 
	else if (dH < -0.25*cells[i][j].ice) dH = -0.25*cells[i][j].ice;
      }
      else {
	if (cells[i-1][j].ice <= 0.0) dH = 0.0;
	else if (dH > 0.25*cells[i-1][j].ice) dH = 0.25*cells[i-1][j].ice;
      }
      vp[i][j].dH = dH;
    }
  }

#pragma omp for schedule(static)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].ice = cells[i][j].oldice + hp[i][j].dH - hp[i][j+1].dH + vp[i][j].dH - vp[i+1][j].dH;
    }
  }
    
  /*compute mean ice thickness and update topice*/
#pragma omp for schedule(static) reduction(+:meanice)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].topice = cells[i][j].bed + cells[i][j].ice;
      meanice += cells[i][j].ice;
    }
  }
#pragma omp single
  {

    meanice /= (double)((*mesh).nc);
    (*mesh).meanice = meanice;
    
  }

  }/*pragma*/

}/*get_change*/

 
void get_change_MUSCL(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,double dt)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double dH,Hminus,Hplus,Hmax,Hmin,rminus,rplus,phi,Hice,a,b;
  double meanice = 0.0; 
  double ss = 1.0e-3;

#pragma omp parallel shared(cells,vp,hp,mesh) private(i,j,dH,Hminus,Hplus,Hmax,Hmin,rminus,rplus,phi,a,b,Hice) firstprivate(nx,ny,dx,dy,dt,ss)
  {
    
    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].oldice = cells[i][j].ice;
      }
    }


  /*loop h-points*/
#pragma omp for schedule(static)
  for (i=0;i<ny;i++) {
    for (j=1;j<nx;j++) {
  
      if (j == 1) Hminus = cells[i][j-1].ice;  
      else {
	rminus = (cells[i][j-1].ice-cells[i][j-2].ice)/(cells[i][j].ice-cells[i][j-1].ice+ss);

	/*minmod limiter*/
	if (rminus > 1.0) rminus = 1.0; if (rminus < 0.0) rminus = 0.0; 

	/*suberbee limiter*/
	/*a = 2.0*rminus; if (a > 1.0) a = 1.0;
	b = rminus; if (b > 2.0) b = 2.0;
	if (a >= b) rminus = a;
	else rminus = b;
	if (rminus < 0.0) rminus = 0.0; */

	Hminus = cells[i][j-1].ice+0.5*rminus*(cells[i][j].ice-cells[i][j-1].ice); 
      }
 
      if (j == (nx-1)) Hplus = cells[i][j].ice; 
      else {
	rplus = (cells[i][j].ice-cells[i][j-1].ice)/(cells[i][j+1].ice-cells[i][j].ice+ss);
	

	/*minmod limiter*/
	if (rplus > 1.0) rplus = 1.0; if (rplus < 0.0) rplus = 0.0;

	/*suberbee limiter*/
	/*a = 2.0*rplus; if (a > 1.0) a = 1.0;
	b = rplus; if (b > 2.0) b = 2.0;
	if (a >= b) rplus = a;
	else rplus = b;
	if (rplus < 0.0) rplus = 0.0; */

	Hplus = cells[i][j].ice-0.5*rplus*(cells[i][j+1].ice-cells[i][j].ice); 

      }
 
      if (Hplus > Hminus) {
	Hmax = Hplus;
	Hmin = Hminus;
      }
      else {
	Hmax = Hminus;
	Hmin = Hplus;
      }
 
      if ((hp[i][j].vx >= 0.0)&&(Hminus <= Hplus)) Hice = Hmin;
      else if ((hp[i][j].vx >= 0.0)&&(Hminus > Hplus)) Hice = Hmax;
      else if ((hp[i][j].vx < 0.0)&&(Hminus <= Hplus)) Hice = Hmax;
      else if ((hp[i][j].vx < 0.0)&&(Hminus > Hplus)) Hice = Hmin;

      dH = hp[i][j].vx*Hice*dt/dx;
      hp[i][j].dH = dH; 
      /*hp[i][j].ice = Hice;*/
    }
  }

  /*loop v-points*/
#pragma omp for schedule(static)
  for (i=1;i<ny;i++) {
    for (j=0;j<nx;j++) { 

      if (i == 1) Hminus = cells[i-1][j].ice; 
      else { 
	rminus = (cells[i-1][j].ice-cells[i-2][j].ice)/(cells[i][j].ice-cells[i-1][j].ice+ss);

	/*minmod limiter*/
	if (rminus > 1.0) rminus = 1.0; if (rminus < 0.0) rminus = 0.0;

	/*suberbee limiter*/
	/*a = 2.0*rminus; if (a > 1.0) a = 1.0;
	b = rminus; if (b > 2.0) b = 2.0;
	if (a >= b) rminus = a;
	else rminus = b;
	if (rminus < 0.0) rminus = 0.0; */

	Hminus = cells[i-1][j].ice+0.5*rminus*(cells[i][j].ice-cells[i-1][j].ice); 
      }
 
      if (i == (ny-1)) Hplus = cells[i][j].ice; 
      else {
	rplus = (cells[i][j].ice-cells[i-1][j].ice)/(cells[i+1][j].ice-cells[i][j].ice+ss);

	/*minmod limiter*/
	if (rplus > 1.0) rplus = 1.0; if (rplus < 0.0) rplus = 0.0;

	/*suberbee limiter*/
	/*a = 2.0*rplus; if (a > 1.0) a = 1.0;
	b = rplus; if (b > 2.0) b = 2.0;
	if (a >= b) rplus = a;
	else rplus = b;
	if (rplus < 0.0) rplus = 0.0; */

	Hplus = cells[i][j].ice-0.5*rplus*(cells[i+1][j].ice-cells[i][j].ice); 
      }

      if (Hplus > Hminus) {
	Hmax = Hplus;
	Hmin = Hminus;
      }
      else {
	Hmax = Hminus;
	Hmin = Hplus;
      }
 
      if ((vp[i][j].vy >= 0.0)&&(Hminus <= Hplus)) Hice = Hmin;
      else if ((vp[i][j].vy >= 0.0)&&(Hminus > Hplus)) Hice = Hmax;
      else if ((vp[i][j].vy < 0.0)&&(Hminus <= Hplus)) Hice = Hmax;
      else if ((vp[i][j].vy < 0.0)&&(Hminus > Hplus)) Hice = Hmin;

      dH = vp[i][j].vy*Hice*dt/dy;
      vp[i][j].dH = dH; 
      /*vp[i][j].ice = Hice;*/
    }
  }

#pragma omp for schedule(static)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].ice = cells[i][j].oldice + hp[i][j].dH - hp[i][j+1].dH + vp[i][j].dH - vp[i+1][j].dH;
    }
  }
    
  /*compute mean ice thickness and update topice*/
#pragma omp for schedule(static) reduction(+:meanice)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].topice = cells[i][j].bed + cells[i][j].ice;
      meanice += cells[i][j].ice;
    }
  }
#pragma omp single
  {

    meanice /= (double)((*mesh).nc);
    (*mesh).meanice = meanice;
    
  }

  }/*pragma*/

}/*get_change_MUSCL*/


void glacial_erosion(celltype **cells,hptype **hp,vptype **vp,iproptype iprop,meshtype *mesh,double dt)
{

  int i,j;
  double ub,slope,pe,qtot=0.0,atot=0.0,meanqrate=0.0,meanarate=0.0;
  double Kq = iprop.Kq;
  double Ka = iprop.Ka;
  double ap = iprop.ap;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int doglacialsedi = (*mesh).doglacialsedi;
  int doparticles = (*mesh).doparticles;
  double efac;
  double minefac = iprop.minefac;
  double sedifac = iprop.sedifac;
  double qa = 1.0;
  double aa = 0.1;
  double maxerate = 1.0e-1;

  
#pragma omp parallel shared(cells,mesh,meanqrate,meanarate) private(i,j,ub,slope,pe,efac) firstprivate(nx,ny,Ka,Kq,dt,qtot,atot,doglacialsedi,doparticles,qa,aa,minefac,sedifac,maxerate)
  {

    /*zero erosion rates*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].quarrying_rate = 0.0;
	cells[i][j].abrasion_rate = 0.0;
      }
    }


    /*loop cells*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if ((cells[i][j].sliding > 0.01)&&(cells[i][j].margin < 0)) {
	  ub = cells[i][j].sliding;
	  slope = -(cells[i][j].dbdx*cells[i][j].vx_b+cells[i][j].dbdy*cells[i][j].vy_b)/ub;
	  pe = 0.01*(cells[i][j].tn-0.25*cells[i][j].Pw);
	  if (pe < 0.0) pe = 0.0; if (pe > 3.0) pe = 3.0;
	  if (slope > -0.55) {
	    if (slope > 1.0) slope = 1.0;
	    /*qtot = pow(pe,3.0)*pow(ub,1.0)*pow(slope+0.55,2.0);*/
	    qtot = pe*ub*(slope+0.55);
	  }
	  else {
	    qtot = 0.0;
	  }
	  qtot *= cells[i][j].Kq*sqrt(1.0+pow(cells[i][j].bslope,2.0));

	  /*
	  if (doglacialsedi > 0) {
	    efac = cells[i][j].cbs*cells[i][j].angular;
	    if (efac < minefac) efac = minefac;  
	    atot = Ka*efac*cells[i][j].ts*ub*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	  }
	  else {
	    atot = Ka*pow(ub,ap)*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	    }*/


	  if (doparticles > 0) {
	    efac = cells[i][j].afac;
	    if (efac < minefac) efac = minefac;   
	    /*if (efac > 10.0) efac = 10.0;*/
	    /*atot = Ka*pow(ub,ap)*sqrt(1.0+pow(cells[i][j].bslope,2.0));*/
	    /*atot = cells[i][j].Ka*efac*cells[i][j].ts*ub*sqrt(1.0+pow(cells[i][j].bslope,2.0));*/
	    atot = cells[i][j].Ka*efac*ub*ub*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	   }
	  else {
	    atot = cells[i][j].Ka*cells[i][j].ts*pow(ub,ap)*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	    }


	  if (atot > maxerate) atot = maxerate;
	  if (qtot > maxerate) qtot = maxerate;


	  /*erode cell*/
	  cells[i][j].quarrying_rate = qtot;
	  cells[i][j].abrasion_rate = atot;
	  cells[i][j].quarrying += qtot*dt;
	  cells[i][j].abrasion += atot*dt; 
	  cells[i][j].bed -= (qtot+atot)*dt; 

	  /*allow particles to transport sediment*/
	  if (doparticles > 0) {
	    /*cells[i][j].sedi += (qtot+atot)*dt; */
	    cells[i][j].sedi += qtot*dt; /*only quarried material stays*/ 
	  }
	  

	  /*add sediment*/
	  if (doglacialsedi > 0) {
	    cells[i][j].angular = cells[i][j].angular*cells[i][j].Vs[19]+sedifac*(qtot*dt*qa+atot*dt*aa);
	    cells[i][j].Vs[19] += sedifac*(qtot+atot)*dt;
	    cells[i][j].angular /= (cells[i][j].Vs[19]+1e-12);
	  }
	  
	}/*if sliding*/
      }/*j*/
    }/*i*/
    
    /*compute mean ratea*/
#pragma omp for schedule(static) reduction(+:meanqrate) reduction(+:meanarate)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	meanqrate += cells[i][j].quarrying_rate;
	meanarate += cells[i][j].abrasion_rate;
      }
    }
#pragma omp single
    {
      meanqrate /= (double)((*mesh).nc);
      meanarate /= (double)((*mesh).nc);
      (*mesh).mean_quarrying_rate = meanqrate;
      (*mesh).mean_abrasion_rate = meanarate;
    }


  }/*pragma*/
}

void get_margin_old(celltype **cells,meshtype *mesh)
{

  int i,j;
  double minice = 10.0;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  long nice = 0;

  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {

      if (cells[i][j].ice < minice) cells[i][j].margin = 0;
      else {

	if (cells[i][j-1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i+1][j].ice < minice) cells[i][j].margin = 1;
	else if (cells[i][j+1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i-1][j].ice < minice) cells[i][j].margin = 1;
	else {
	  cells[i][j].margin = -1; 
	  nice += 1;
	}

      }/*else*/

    }/*j*/
  }/*i*/
  (*mesh).nice = nice;
}

void get_margin(celltype **cells,meshtype *mesh)
{

  int i,j;
  double minice = 10.0;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  long nice = 0;

  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {

      if (cells[i][j].ice < minice) {
       
	cells[i][j].margin = 0;

	if ((cells[i][j].topsedi < cells[i-1][j].topsedi)&&(cells[i][j].topsedi < cells[i+1][j].topsedi)&&(cells[i][j].topsedi < cells[i][j-1].topsedi)&&(cells[i][j].topsedi < cells[i][j+1].topsedi))
	  cells[i][j].pit = 1;
	else cells[i][j].pit = 0;

      }
      else {

	if (cells[i][j-1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i+1][j].ice < minice) cells[i][j].margin = 1;
	else if (cells[i][j+1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i-1][j].ice < minice) cells[i][j].margin = 1;
	else {
	  cells[i][j].margin = -1; 
	  nice += 1;
	}

      }/*else*/

    }/*j*/
  }/*i*/

  /*boundaries*/
  for (i=0;i<ny;i++) {
    if (cells[i][0].ice > minice) cells[i][0].margin = -1;
    else cells[i][0].margin = 0;
    if (cells[i][nx-1].ice > minice) cells[i][nx-1].margin = -1;
    else cells[i][nx-1].margin = 0;
  }

  for (j=0;j<nx;j++) {
    if (cells[0][j].ice > minice) cells[0][j].margin = -1;
    else cells[0][j].margin = 0;
    if (cells[ny-1][j].ice > minice) cells[ny-1][j].margin = -1;
    else cells[ny-1][j].margin = 0;
  }


  (*mesh).nice = nice;
}


void compute_leefactor(celltype **cells,meshtype mesh)
{

  int i,j,k;
  int w;
  int wd = 20;
  double ss = 0.1;
  double s,maxs1,maxs2,maxs3,maxs4;

  int nx = mesh.nx;
  int ny = mesh.ny;
  double dx = mesh.dx;
  double dy = mesh.dy;

  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      if (cells[i][j].ice > 100.0) cells[i][j].lee = 3.0;
      else {
      
	maxs1 = 0.0;
	if (j > wd) w = wd; else w = j;
	for (k=1;k<=w;k++) {
	  s = (cells[i][j-k].bed - cells[i][j].topice)/(k*dx);
	  if (s > maxs1) maxs1 = s;
	}
	
	maxs2 = 0.0;
	if (j < nx-wd-1) w = wd; else w = nx-j-1;
	for (k=1;k<=w;k++) {
	  s = (cells[i][j+k].bed - cells[i][j].topice)/(k*dx);
	  if (s > maxs2) maxs2 = s;
	}
	

	maxs3 = 0.0;
	if (i > wd) w = wd; else w = i;
	for (k=1;k<=w;k++) {
	  s = (cells[i-k][j].bed - cells[i][j].topice)/(k*dy);
	  if (s > maxs3) maxs3 = s;
	}
	
	maxs4 = 0.0;
	if (i < ny-wd-1) w = wd; else w = ny-i-1;
	for (k=1;k<=w;k++) {
	  s = (cells[i+k][j].bed - cells[i][j].topice)/(k*dy);
	  if (s > maxs4) maxs4 = s;
	}
	

	maxs1 = maxs1/ss; if (maxs1 > 1.0) maxs1 = 1.0;
	maxs2 = maxs2/ss; if (maxs2 > 1.0) maxs2 = 1.0;
	maxs3 = maxs3/ss; if (maxs3 > 1.0) maxs3 = 1.0;
	maxs4 = maxs4/ss; if (maxs4 > 1.0) maxs4 = 1.0;
	
	/*
	  maxs1 = 1.0 - maxs1;
	  maxs2 = 1.0 - maxs2;
	  maxs3 = 1.0 - maxs3;
	  maxs4 = 1.0 - maxs4;
	*/
	
	cells[i][j].lee = maxs1+maxs2+maxs3+maxs4-1.0;
	if (cells[i][j].lee < 0.0) cells[i][j].lee = 0.0;
	/*cells[i][j].lee = maxs1*maxs2*maxs3*maxs4;*/
      }

    }
  }


}

void glacial_sediment_transport(celltype **cells,hptype **hp,vptype **vp,meshtype mesh,iproptype iprop,double dt)
{

  /*array Vs stores sediment thickness in a 3D grid*/
  /*Vs[0] is uppermost layer*/

  int i,j,k,t;
  double dz,zz,ice;
  double ms,mb,dc,cs,ds;
  double vx,vy,vz;

  int ny = mesh.ny;
  int nx = mesh.nx;
  double dx = mesh.dx;
  double dy = mesh.dy;
  int dosediment = mesh.dosediment;

  int Nz = 20;
  double minice = 20.0;
  double La = 3.0e3;
  double ksg = iprop.ksg*iprop.rho*mesh.gravity;

  int Nzt = 10;
  double dtz = (double)(dt/Nzt);

  /*loop cells for vertical advection*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      ice = cells[i][j].ice;

      if (ice > minice) {
	dz = (double)(ice/Nz);
	for (t=0;t<Nzt;t++) {
	  for (k=1;k<Nz;k++) {
	    zz = (double)k*dz;
	    ms = cells[i][j].Ms;
	    mb = -cells[i][j].Mb;
	    vz = (ice-zz)/ice*ms+(zz/ice)*mb; /*vz positive downwards*/
	    if (vz < 0.0) {
	      dc = dtz*vz*cells[i][j].Vs[k]/dz; /*dc<0*/
	      if (dc < -cells[i][j].Vs[k]) dc = -cells[i][j].Vs[k];
	      cells[i][j].Vs[k] += dc;
	      cells[i][j].Vs[k-1] -= dc;
	    }
	    else if (vz > 0.0) {
	      dc = dtz*vz*cells[i][j].Vs[k-1]/dz; /*dc>0*/
	      if (dc > cells[i][j].Vs[k-1]) dc = cells[i][j].Vs[k-1];
	      cells[i][j].Vs[k-1] -= dc;
	      cells[i][j].Vs[k] += dc;
	    }
	  }/*k*/
	}/*t*/
      }/*if*/
    }/*j*/
  }/*i*/

  /*loop hpoints for horizontal advection*/
  for (i=0;i<ny;i++) {
    for (j=1;j<nx;j++) {
      if (hp[i][j].ice > minice) {
	vx = hp[i][j].vx_b + hp[i][j].vx_d; /*uniform velocity*/
	if (vx < 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vx*cells[i][j].Vs[k]/dx;/*dc<0*/
	    if (dc < -cells[i][j].Vs[k]) dc = -cells[i][j].Vs[k];
	    cells[i][j].Vs[k] += dc;
	    cells[i][j-1].Vs[k] -= dc;
	    if (k == (Nz-1)) {
	      cells[i][j-1].angular = (cells[i][j-1].angular*(cells[i][j-1].Vs[k]+dc)-cells[i][j].angular*dc*exp(-fabs(dt*vx)/La))/(cells[i][j-1].Vs[k]+1.0e-6);
	    }
	  }/*k*/
	}/*if*/
	else if (vx > 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vx*cells[i][j-1].Vs[k]/dx;/*dc>0*/
	    if (dc > cells[i][j-1].Vs[k]) dc = cells[i][j-1].Vs[k];
	    cells[i][j-1].Vs[k] -= dc;
	    cells[i][j].Vs[k] += dc;
	    if (k == (Nz-1)) {
	      cells[i][j].angular = (cells[i][j].angular*(cells[i][j].Vs[k]-dc)+cells[i][j-1].angular*dc*exp(-fabs(dt*vx)/La))/(cells[i][j].Vs[k]+1.0e-6);
	    }
	  }/*k*/
	}/*if*/
      }/*if*/
    }/*j*/
  }/*i*/

  /*loop vpoints for horizontal advection*/
  for (i=1;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if (vp[i][j].ice > minice) {
	vy = vp[i][j].vy_b + vp[i][j].vy_d; /*uniform velocity*/
	if (vy < 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vy*cells[i][j].Vs[k]/dy; /*dc<0*/
	    if (dc < -cells[i][j].Vs[k]) dc = -cells[i][j].Vs[k];
	    cells[i][j].Vs[k] += dc;
	    cells[i-1][j].Vs[k] -= dc;
	    if (k == (Nz-1)) {
	      cells[i-1][j].angular = (cells[i-1][j].angular*(cells[i-1][j].Vs[k]+dc)-cells[i][j].angular*dc*exp(-fabs(dt*vy)/La))/(cells[i-1][j].Vs[k]+1.0e-6);
	    }
	  }/*k*/
	}/*if*/
	else if (vy > 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vy*cells[i-1][j].Vs[k]/dy; /*dc>0*/
	    if (dc > cells[i-1][j].Vs[k]) dc = cells[i-1][j].Vs[k];
	    cells[i-1][j].Vs[k] -= dc;
	    cells[i][j].Vs[k] += dc;
	    if (k == (Nz-1)) {
	      cells[i][j].angular = (cells[i][j].angular*(cells[i][j].Vs[k]-dc)+cells[i-1][j].angular*dc*exp(-fabs(dt*vy)/La))/(cells[i][j].Vs[k]+1.0e-6);
	    }
	  }/*k*/
	}/*if*/
      }/*if*/
    }/*j*/
  }/*i*/
   
  /*entrain and deposit sediment at base of ice*/ 
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if ((dosediment > 0)&&(cells[i][j].ice > minice)) {
	dz = (double)(cells[i][j].ice/Nz);
	cells[i][j].vsb = ksg*cells[i][j].te/(cells[i][j].Vs[Nz-1]+0.1)-cells[i][j].Vs[Nz-1]*cells[i][j].Mb/dz;
	ds = cells[i][j].vsb*dt; 
	if (ds > cells[i][j].sedi) ds = cells[i][j].sedi;
	if (ds < -cells[i][j].Vs[Nz-1]) ds = -cells[i][j].Vs[Nz-1];
	cells[i][j].Vs[Nz-1] += ds;
	cells[i][j].sedi -= ds;
	cells[i][j].cbs = cells[i][j].Vs[Nz-1]/dz; /*basal concentration*/
	if (cells[i][j].cbs > 0.2) cells[i][j].cbs = 0.2;
      }
      else {
	cells[i][j].vsb = 0.0;
	cells[i][j].cbs = 0.0;
      }
    }
  }

  /*transfer sediment to sedi array at ice margin*/
  if (dosediment > 0) {
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) { 
	if ((cells[i][j].margin == 1)&&(cells[i][j].ice > minice)) { 
	  dz = (double)(cells[i][j].ice/Nz);
	  for (k=0;k<Nz;k++) {
	    cs = cells[i][j].Vs[k]/dz;
	    if (cs > 0.5) {
	      cells[i][j].sedi += (cs-0.5)*dz;
	      cells[i][j].Vs[k] = 0.5*dz;
	    }
	  }/*k*/ 
	}/*if*/ 
	else if (cells[i][j].ice < minice) {
	  for (k=0;k<Nz;k++) {
	    cells[i][j].sedi += cells[i][j].Vs[k];
	    cells[i][j].Vs[k] = 0.0;
	  } 
	}
      }/*j*/ 
    }/*i*/

  }



}
