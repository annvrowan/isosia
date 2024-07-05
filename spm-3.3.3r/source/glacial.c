
void isosia(celltype **cells,hptype **hp,vptype **vp,cornertype **cp,meshtype *mesh,iproptype iprop)
{

  int i,j;
  long stressit,veloitt;
  double sres,sres_old,sdiff,stot,vres,vres_old,vdiff,vtot;
  double tx,ty,ts,vb,beta,vx,vy,Res,lmaxRes,maxRes = 1.0;
  double dtfac,lmaxdt,maxdt;
  double maxspeed=0.0,lmaxspeed=0.0;
  double minice = 5.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int periodic = (*mesh).periodic;
  int dosliding = (*mesh).dosliding;

  int smode = (*mesh).slidingmode;
  int coldbased = (*mesh).coldbased;
  int gmode = 1;

  double gamma = iprop.gamma;
  double gamma0 = iprop.gamma0;
  double ifac = iprop.ifac;
  double sfac = iprop.sfac;
  double vbfac = iprop.vbfac;
  int maxitt_v = iprop.maxitt_v;
  int maxitt_s = iprop.maxitt_s;
  double Cs = iprop.Cs;
  double maxdef = iprop.maxdeformation;

  double C = iprop.C;
  double L0 = iprop.L0;
  double maxsliding = iprop.maxsliding;
  double minbeta = iprop.minbeta;

  long ndiv = 0;

  /*initialize*/
  maxspeed = 0.1;
  lmaxdt = (*mesh).maxdt;
  maxdt = (*mesh).maxdt;



#pragma omp parallel shared(cells,vp,hp,cp,mesh,stressit,sres,sres_old,sdiff,stot,veloitt,vres,vres_old,vdiff,vtot,maxspeed,maxRes,maxdt,ndiv) private(i,j,tx,ty,ts,vb,beta,vx,vy,Res,lmaxRes,dtfac) firstprivate(nx,ny,dx,dy,periodic,lmaxspeed,gamma,gamma0,Cs,ifac,sfac,vbfac,maxitt_v,maxitt_s,C,L0,minbeta,maxsliding,dosliding,smode,gmode,maxdef,minice,lmaxdt)
  {


#pragma omp single
    {
      veloitt = 0;
      vres = 1.0;
    }

#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  cells[i][j].te2i = cells[i][j].te2;
	}
      }
      
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) {
	  hp[i][j].vx_di = hp[i][j].vx_d;
	  hp[i][j].vx_bi = hp[i][j].vx_b;

	}
      }
#pragma omp for schedule(static)
      for (i=1;i<ny;i++) {
	for (j=0;j<nx;j++) {
	    vp[i][j].vy_di = vp[i][j].vy_d;
	    vp[i][j].vy_bi = vp[i][j].vy_b;
	}
      }


    /*itterate velocities*/
      while ((vres > 1.0e-3)&&(veloitt < maxitt_v)) {
      /*while ((maxRes > 1.0e-2)&&(veloitt < maxitt_v)) {*/
      
      /******* Strain rates ********/
      /*shear strain rates*/
#pragma omp for schedule(static)
	for (i=1;i<ny;i++) {
	  for (j=1;j<nx;j++) {
	    cp[i][j].exy = 0.5*((hp[i][j].vx-hp[i-1][j].vx)/dy + (vp[i][j].vy-vp[i][j-1].vy)/dx);
	  }
	}

      if (periodic > 0) {

	/*shear strain boundary conditions*/
#pragma omp single
	{
	  cp[0][0].exy = 0.5*((hp[0][0].vx-hp[ny-1][0].vx)/dy+(vp[0][0].vy-vp[0][nx-1].vy)/dx);
	  cp[ny][0].exy = 0.5*((hp[0][0].vx-hp[ny-1][0].vx)/dy+(vp[ny][0].vy-vp[ny][nx-1].vy)/dx);
	  cp[ny][nx].exy = 0.5*((hp[0][nx].vx-hp[ny-1][nx].vx)/dy+(vp[ny][0].vy-vp[ny][nx-1].vy)/dx);
	  cp[0][nx].exy = 0.5*((hp[0][nx].vx-hp[nx-1][nx].vx)/dy+(vp[0][0].vy-vp[0][nx-1].vy)/dx);
	}
#pragma omp for schedule(static)
	for (i=1;i<ny;i++) {
	  cp[i][0].exy = 0.5*((hp[i][0].vx-hp[i-1][0].vx)/dy+(vp[i][0].vy-vp[i][nx-1].vy)/dx);
	  cp[i][nx].exy = cp[i][0].exy;
	}
#pragma omp for schedule(static)
	for (j=1;j<nx;j++) {
	  cp[0][j].exy = 0.5*((hp[0][j].vx-hp[ny-1][j].vx)/dy+(vp[0][j].vy-vp[0][j-1].vy)/dx);
	  cp[ny][j].exy = cp[0][j].exy;
	}
      }      


      /*strain rates in cells*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  if ((cells[i][j].ice > minice)&&(gmode > 0)) {
	    cells[i][j].exx = (hp[i][j+1].vx-hp[i][j].vx)/dx-(cells[i][j].dhdx*((hp[i][j].vx_s-hp[i][j].vx)+(hp[i][j+1].vx_s-hp[i][j+1].vx))/2.0+cells[i][j].dbdx*((hp[i][j].vx-hp[i][j].vx_b)+(hp[i][j+1].vx-hp[i][j+1].vx_b))/2.0)/cells[i][j].ice;
	    cells[i][j].eyy = (vp[i+1][j].vy-vp[i][j].vy)/dy-(cells[i][j].dhdy*((vp[i][j].vy_s-vp[i][j].vy)+(vp[i+1][j].vy_s-vp[i+1][j].vy))/2.0+cells[i][j].dbdy*((vp[i][j].vy-vp[i][j].vy_b)+(vp[i+1][j].vy-vp[i+1][j].vy_b))/2.0)/cells[i][j].ice;
	    cells[i][j].exy = 0.25*(cp[i][j].exy+cp[i+1][j].exy+cp[i+1][j+1].exy+cp[i][j+1].exy)-(cells[i][j].dhdy*((hp[i][j].vx_s-hp[i][j].vx)+(hp[i][j+1].vx_s-hp[i][j+1].vx))/2.0+cells[i][j].dhdx*((vp[i][j].vy_s-vp[i][j].vy)+(vp[i+1][j].vy_s-vp[i+1][j].vy))/2.0+cells[i][j].dbdy*((hp[i][j].vx-hp[i][j].vx_b)+(hp[i][j+1].vx-hp[i][j+1].vx_b))/2.0+cells[i][j].dbdx*((vp[i][j].vy-vp[i][j].vy_b)+(vp[i+1][j].vy-vp[i+1][j].vy_b))/2.0)/(2.0*cells[i][j].ice);
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
	      cells[i][j].sxx = cells[i][j].exx/(gamma*cells[i][j].te2+gamma0);
	      cells[i][j].syy = cells[i][j].eyy/(gamma*cells[i][j].te2+gamma0);
	      cells[i][j].sxy = cells[i][j].exy/(gamma*cells[i][j].te2+gamma0);
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
	  hp[i][0].dsxxdx = hp[i][1].dsxxdx;
	  hp[i][0].dsyydx = hp[i][1].dsyydx;
	  hp[i][0].dsxydx = hp[i][1].dsxydx;
	  hp[i][nx].dsxxdx = hp[i][nx-1].dsxxdx;
	  hp[i][nx].dsyydx = hp[i][nx-1].dsyydx;
	  hp[i][nx].dsxydx = hp[i][nx-1].dsxydx;
	}
	
	/*compute stress gradients - y dir*/
#pragma omp for schedule(static)
	for (i=1;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    vp[i][j].dsxxdy = (cells[i][j].sxx-cells[i-1][j].sxx)/dy;
	    vp[i][j].dsyydy = (cells[i][j].syy-cells[i-1][j].syy)/dy;
	    vp[i][j].dsxydy = (cells[i][j].sxy-cells[i-1][j].sxy)/dy;
	  }
	}
	
	if (periodic > 0) {
	
	  /*boundary conditions - stress gradients*/
#pragma omp for schedule(static)
	  for (i=0;i<ny;i++) {
	    hp[i][0].dsxxdx = (cells[i][0].sxx-cells[i][nx-1].sxx)/dx;
	    hp[i][0].dsyydx = (cells[i][0].syy-cells[i][nx-1].syy)/dx;
	    hp[i][0].dsxydx = (cells[i][0].sxy-cells[i][nx-1].sxy)/dx;
	    hp[i][nx].dsxxdx = hp[i][0].dsxxdx;
	    hp[i][nx].dsyydx = hp[i][0].dsyydx;
	    hp[i][nx].dsxydx = hp[i][0].dsxydx;
	  }
#pragma omp for schedule(static)
	  for (j=0;j<nx;j++) {
	    vp[0][j].dsxxdy = (cells[0][j].sxx-cells[ny-1][j].sxx)/dy;
	    vp[0][j].dsyydy = (cells[0][j].syy-cells[ny-1][j].syy)/dy;
	    vp[0][j].dsxydy = (cells[0][j].sxy-cells[ny-1][j].sxy)/dy;
	    vp[ny][j].dsxxdy = vp[0][j].dsxxdy;
	    vp[ny][j].dsyydy = vp[0][j].dsyydy;
	    vp[ny][j].dsxydy = vp[0][j].dsxydy;
	  }
	}	


	/*transfer stress gradients compute coeeficients and update stress*/
#pragma omp for schedule(static) reduction(+:sdiff) reduction(+:stot)
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    cells[i][j].dsxxdx = 0.5*(hp[i][j].dsxxdx+hp[i][j+1].dsxxdx);
	    cells[i][j].dsyydx = 0.5*(hp[i][j].dsyydx+hp[i][j+1].dsyydx);
	    cells[i][j].dsxydx = 0.5*(hp[i][j].dsxydx+hp[i][j+1].dsxydx);
	    cells[i][j].dsxxdy = 0.5*(vp[i][j].dsxxdy+vp[i+1][j].dsxxdy);
	    cells[i][j].dsyydy = 0.5*(vp[i][j].dsyydy+vp[i+1][j].dsyydy);
	    cells[i][j].dsxydy = 0.5*(vp[i][j].dsxydy+vp[i+1][j].dsxydy);
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

	/*	if ((sres > sres_old)&&(sfac > 1.0e-4)) {

#pragma omp for schedule(static)
	  for (i=0;i<ny;i++) {
	    for (j=0;j<nx;j++) {
	      cells[i][j].te2 = cells[i][j].te2i;
	    }
	  }

#pragma omp single 
	  {
	    sfac /= 2.0;
	    stressit = 0;
	    printf("   stressit = %ld, sres = %4.4e, sfac = %4.4e\n",stressit,sres,sfac);
	  }

	}*/
	
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
	    hp[i][j].cx[0] = (cells[i][j-1].sxx+cells[i][j].sxx+.5*(cells[i][j-1].syy+cells[i][j].syy))*hp[i][j].dhdx+.5*(cells[i][j-1].sxy+cells[i][j].sxy)*hp[i][j].dhdy;
	    hp[i][j].cx[1] = -(1.0+hp[i][j].alpha2)*hp[i][j].dhdx+2.0*hp[i][j].dsxxdx+hp[i][j].dsyydx+.5*(cells[i][j-1].dsxydy+cells[i][j].dsxydy);
	    hp[i][j].cy[0] = (.5*(cells[i][j-1].sxx+cells[i][j].sxx)+cells[i][j-1].syy+cells[i][j].syy)*hp[i][j].dhdy+.5*(cells[i][j-1].sxy+cells[i][j].sxy)*hp[i][j].dhdx;
	    hp[i][j].cy[1] = -(1.0+hp[i][j].alpha2)*hp[i][j].dhdy+.5*(cells[i][j-1].dsxxdy+cells[i][j].dsxxdy)+cells[i][j-1].dsyydy+cells[i][j].dsyydy+hp[i][j].dsxydx;
	    hp[i][j].kc[0] = hp[i][j].cx[0]*hp[i][j].cx[0]+hp[i][j].cy[0]*hp[i][j].cy[0]+.5*(cells[i][j-1].td2+cells[i][j].td2);
	    hp[i][j].kc[1] = 2.0*(hp[i][j].cx[0]*hp[i][j].cx[1]+hp[i][j].cy[0]*hp[i][j].cy[1]);
	    hp[i][j].kc[2] = hp[i][j].cx[1]*hp[i][j].cx[1]+hp[i][j].cy[1]*hp[i][j].cy[1];
	    hp[i][j].wx[0] = hp[i][j].cx[0]*hp[i][j].kc[0];/*-((hp[i][j+1].vz_b-hp[i][j-1].vz_b)/dx-hp[i][j].dbdx*(cells[i][j-1].ezz+cells[i][j].ezz)/2.0)/(2.0*gamma);*/
	    hp[i][j].wx[1] = hp[i][j].cx[0]*hp[i][j].kc[1]+hp[i][j].cx[1]*hp[i][j].kc[0];/*-1.0/(4.0*gamma)*(cells[i][j].ezz-cells[i][j-1].ezz)/dx;*/
	    hp[i][j].wx[2] = hp[i][j].cx[0]*hp[i][j].kc[2]+hp[i][j].cx[1]*hp[i][j].kc[1];
	    hp[i][j].wx[3] = hp[i][j].cx[1]*hp[i][j].kc[2];
	    hp[i][j].vrx = 2.0*gamma*(hp[i][j].wx[0]*hp[i][j].ice/2.0+hp[i][j].wx[1]*pow(hp[i][j].ice,2.0)/3.0+hp[i][j].wx[2]*pow(hp[i][j].ice,3.0)/4.0+hp[i][j].wx[3]*pow(hp[i][j].ice,4.0)/5.0);
	    if (fabs(hp[i][j].vrx) > maxdef) hp[i][j].vrx *= maxdef/fabs(hp[i][j].vrx);
	    hp[i][j].vx_d += ifac*(hp[i][j].vrx - hp[i][j].vx_d);
	    Res = fabs(hp[i][j].vrx-hp[i][j].vx_d)/(fabs(hp[i][j].vx_d)+1.0);
	    if (Res > lmaxRes) lmaxRes = Res;
	    hp[i][j].vresx = fabs(hp[i][j].vrx - hp[i][j].vx_d);
	    vdiff += pow(hp[i][j].vrx-hp[i][j].vx_d,2.0);
	    vtot += pow(hp[i][j].vrx,2.0);
	  }
	  else {
	    hp[i][j].vx_d = 0.0;
	    hp[i][j].vresx = 0.0;;
	  }
	}
	hp[i][0].vx_d = hp[i][1].vx_d;
	hp[i][nx].vx_d = hp[i][nx-1].vx_d;
      }

      
      /*y velocity*/
#pragma omp for schedule(static) reduction(+:vdiff) reduction(+:vtot)
      for (i=1;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  if (vp[i][j].ice > minice) {
	    vp[i][j].cx[0] = (cells[i-1][j].sxx+cells[i][j].sxx+.5*(cells[i-1][j].syy+cells[i][j].syy))*vp[i][j].dhdx+.5*(cells[i-1][j].sxy+cells[i][j].sxy)*vp[i][j].dhdy;
	    vp[i][j].cx[1] = -(1.0+vp[i][j].alpha2)*vp[i][j].dhdx+(cells[i-1][j].dsxxdx+cells[i][j].dsxxdx)+.5*(cells[i-1][j].dsyydx+cells[i][j].dsyydx)+vp[i][j].dsxydy;
	    vp[i][j].cy[0] = (.5*(cells[i-1][j].sxx+cells[i][j].sxx)+cells[i-1][j].syy+cells[i][j].syy)*vp[i][j].dhdy+.5*(cells[i-1][j].sxy+cells[i][j].sxy)*vp[i][j].dhdx;
	    vp[i][j].cy[1] = -(1.0+vp[i][j].alpha2)*vp[i][j].dhdy+vp[i][j].dsxxdy+2.0*vp[i][j].dsyydy+.5*(cells[i-1][j].dsxydx+cells[i][j].dsxydx);
	    vp[i][j].kc[0] = vp[i][j].cx[0]*vp[i][j].cx[0]+vp[i][j].cy[0]*vp[i][j].cy[0]+.5*(cells[i-1][j].td2+cells[i][j].td2);
	    vp[i][j].kc[1] = 2.0*(vp[i][j].cx[0]*vp[i][j].cx[1]+vp[i][j].cy[0]*vp[i][j].cy[1]);
	    vp[i][j].kc[2] = vp[i][j].cx[1]*vp[i][j].cx[1]+vp[i][j].cy[1]*vp[i][j].cy[1];
	    vp[i][j].wy[0] = vp[i][j].cy[0]*vp[i][j].kc[0];/*-((vp[i+1][j].vz_b-vp[i-1][j].vz_b)/dy-vp[i][j].dbdy*(cells[i-1][j].ezz+cells[i][j].ezz)/2.0)/(2.0*gamma);*/
	    vp[i][j].wy[1] = vp[i][j].cy[0]*vp[i][j].kc[1]+vp[i][j].cy[1]*vp[i][j].kc[0];/*-1.0/(4.0*gamma)*(cells[i][j].ezz-cells[i-1][j].ezz)/dy;*/
	    vp[i][j].wy[2] = vp[i][j].cy[0]*vp[i][j].kc[2]+vp[i][j].cy[1]*vp[i][j].kc[1];
	    vp[i][j].wy[3] = vp[i][j].cy[1]*vp[i][j].kc[2];
	    vp[i][j].vry = 2.0*gamma*(vp[i][j].wy[0]*vp[i][j].ice/2.0+vp[i][j].wy[1]*pow(vp[i][j].ice,2.0)/3.0+vp[i][j].wy[2]*pow(vp[i][j].ice,3.0)/4.0+vp[i][j].wy[3]*pow(vp[i][j].ice,4.0)/5.0);
	    if (fabs(vp[i][j].vry) > maxdef) vp[i][j].vry *= maxdef/fabs(vp[i][j].vry);
	    vp[i][j].vy_d += ifac*(vp[i][j].vry - vp[i][j].vy_d);
	    Res = fabs(vp[i][j].vry-vp[i][j].vy_d)/(fabs(vp[i][j].vy_d)+1.0);
	    if (Res > lmaxRes) lmaxRes = Res;
	    vp[i][j].vresy = fabs(vp[i][j].vry - vp[i][j].vy_d);
	    vdiff += pow(vp[i][j].vry-vp[i][j].vy_d,2.0);
	    vtot += pow(vp[i][j].vry,2.0);
	  }
	  else {
	    vp[i][j].vy_d = 0.0;
	    vp[i][j].vresy = 0.0;
	  }
	}
	vp[0][i].vy_d = vp[1][i].vy_d;
	vp[ny][i].vy_d = vp[ny-1][i].vy_d;
      }
#pragma omp critical 
      { 
	if (lmaxRes > maxRes) maxRes = lmaxRes; 
      }


      if (periodic > 0) {
      
	/*boundary conditions*/
#pragma omp for schedule(static) reduction(+:vdiff) reduction(+:vtot)
	for (i=0;i<ny;i++) {
	  hp[i][0].cx[0] = 0.5*(cells[i][0].cx[0]+cells[i][nx-1].cx[0]);
	  hp[i][0].cx[1] = 0.5*(cells[i][0].cx[1]+cells[i][nx-1].cx[1]);
	  hp[i][0].kc[0] = 0.5*(cells[i][0].kc[0]+cells[i][nx-1].kc[0]);
	  hp[i][0].kc[1] = 0.5*(cells[i][0].kc[1]+cells[i][nx-1].kc[1]);
	  hp[i][0].kc[2] = 0.5*(cells[i][0].kc[2]+cells[i][nx-1].kc[2]);
	  hp[i][0].wx[0] = hp[i][0].cx[0]*hp[i][0].kc[0];
	  hp[i][0].wx[1] = hp[i][0].cx[0]*hp[i][0].kc[1]+hp[i][0].cx[1]*hp[i][0].kc[0];
	  hp[i][0].wx[2] = hp[i][0].cx[0]*hp[i][0].kc[2]+hp[i][0].cx[1]*hp[i][0].kc[1];
	  hp[i][0].wx[3] = hp[i][0].cx[1]*hp[i][0].kc[2];
	  hp[i][0].vrx = 2.0*gamma*(hp[i][0].wx[0]*hp[i][0].ice/2.0+hp[i][0].wx[1]*pow(hp[i][0].ice,2.0)/3.0+hp[i][0].wx[2]*pow(hp[i][0].ice,3.0)/4.0+hp[i][0].wx[3]*pow(hp[i][0].ice,4.0)/5.0);
	  hp[i][0].vx_d += ifac*(hp[i][0].vrx - hp[i][0].vx_d);
	  vdiff += pow(hp[i][0].vrx-hp[i][0].vx_d,2.0);
	  vtot += pow(hp[i][0].vrx,2.0);
	  hp[i][nx].vx_d = hp[i][0].vx_d;
	}
#pragma omp for schedule(static) reduction(+:vdiff) reduction(+:vtot)
	for (j=0;j<nx;j++) {
	  vp[0][j].cy[0] = 0.5*(cells[0][j].cy[0]+cells[ny-1][j].cy[0]);
	  vp[0][j].cy[1] = 0.5*(cells[0][j].cy[1]+cells[ny-1][j].cy[1]);
	  vp[0][j].kc[0] = 0.5*(cells[0][j].kc[0]+cells[ny-1][j].kc[0]);
	  vp[0][j].kc[1] = 0.5*(cells[0][j].kc[1]+cells[ny-1][j].kc[1]);
	  vp[0][j].kc[2] = 0.5*(cells[0][j].kc[2]+cells[ny-1][j].kc[2]);
	  vp[0][j].wy[0] = vp[0][j].cy[0]*vp[0][j].kc[0];
	  vp[0][j].wy[1] = vp[0][j].cy[0]*vp[0][j].kc[1]+vp[0][j].cy[1]*vp[0][j].kc[0];
	  vp[0][j].wy[2] = vp[0][j].cy[0]*vp[0][j].kc[2]+vp[0][j].cy[1]*vp[0][j].kc[1];
	  vp[0][j].wy[3] = vp[0][j].cy[1]*vp[0][j].kc[2];
	  vp[0][j].vry = 2.0*gamma*(vp[0][j].wy[0]*vp[0][j].ice/2.0+vp[0][j].wy[1]*pow(vp[0][j].ice,2.0)/3.0+vp[0][j].wy[2]*pow(vp[0][j].ice,3.0)/4.0+vp[0][j].wy[3]*pow(vp[0][j].ice,4.0)/5.0);
	  vp[0][j].vy_d += ifac*(vp[0][j].vry - vp[0][j].vy_d);
	  vdiff += pow(vp[0][j].vry-vp[0][j].vy_d,2.0);
	  vtot += pow(vp[0][j].vry,2.0);
	  vp[ny][j].vy = vp[0][j].vy;
	}
      
      }
      

      /******** basal sliding **********/
      if (dosliding > 0) {

#pragma omp for schedule(static) 
	for (i=0;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    cells[i][j].lb = sqrt(1.0+pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0));
	    cells[i][j].pb = (1.0+pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0))*cells[i][j].ice+cells[i][j].szz;
	    cells[i][j].sxz = cells[i][j].cx[0]+cells[i][j].cx[1]*cells[i][j].ice;
	    cells[i][j].syz = cells[i][j].cy[0]+cells[i][j].cy[1]*cells[i][j].ice;
	    cells[i][j].tn = cells[i][j].pb+(2.0*cells[i][j].dbdx*cells[i][j].sxz+2.0*cells[i][j].dbdy*cells[i][j].syz-pow(cells[i][j].dbdx,2.0)*cells[i][j].sxx-pow(cells[i][j].dbdy,2.0)*cells[i][j].syy-cells[i][j].szz-2.0*cells[i][j].dbdx*cells[i][j].dbdy*cells[i][j].sxy)/pow(cells[i][j].lb,2.0);
	    /*cells[i][j].tn = cells[i][j].ice;*/
	    if (cells[i][j].tn < 0.0) cells[i][j].tn = 0.0;
	    cells[i][j].te = cells[i][j].tn-cells[i][j].Pw; if (cells[i][j].te < 0.0) cells[i][j].te = 0.0;/* if (cells[i][j].te > 300.0) cells[i][j].te = 300.0;*/
	    cells[i][j].tbx = (cells[i][j].dbdx*(cells[i][j].pb-cells[i][j].sxx)-cells[i][j].dbdy*cells[i][j].sxy+cells[i][j].sxz-cells[i][j].dbdx*cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].tby = (-cells[i][j].dbdx*cells[i][j].sxy+cells[i][j].dbdy*(cells[i][j].pb-cells[i][j].syy)+cells[i][j].syz-cells[i][j].dbdy*cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].tbz = (-cells[i][j].dbdx*cells[i][j].sxz-cells[i][j].dbdy*cells[i][j].syz-cells[i][j].pb+cells[i][j].szz+cells[i][j].tn)/cells[i][j].lb;
	    cells[i][j].ts = sqrt(pow(cells[i][j].tbx,2.0)+pow(cells[i][j].tby,2.0)+pow(cells[i][j].tbz,2.0));
	    if (smode == 0) {
	      if (cells[i][j].ice > 10.0) cells[i][j].vb = Cs*pow(cells[i][j].ts,2.0);
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (Cs*cells[i][j].ice/10.0)*pow(cells[i][j].ts,2.0);
	      else cells[i][j].vb = 0.0;
	    }
	    else if (smode == 1) {
	      if (cells[i][j].ice > 10.0) cells[i][j].vb = Cs*pow(cells[i][j].ts,2.0)/(cells[i][j].te+1.0);
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (Cs*cells[i][j].ice/10.0)*pow(cells[i][j].ts,2.0)/(cells[i][j].te+1.0);
	      else cells[i][j].vb = 0.0;
	    }
	    else {
	      beta = pow(C,3.0)-pow(cells[i][j].ts/(cells[i][j].te+1.0e-3),3.0);
	      if (beta < minbeta) beta = minbeta;
	      if (beta <= L0*pow(cells[i][j].ts,3.0)/maxsliding) cells[i][j].vb = maxsliding;
	      else if (cells[i][j].ice > 10.0) cells[i][j].vb = L0*pow(cells[i][j].ts,3.0)/beta;
	      else if (cells[i][j].ice > 1.0) cells[i][j].vb = (L0*cells[i][j].ice/10.0)*pow(cells[i][j].ts,3.0)/beta;
	      else cells[i][j].vb = 0.0;	    
	    }
	    if (coldbased == 1) cells[i][j].vb *= cells[i][j].sfac;
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
	  hp[i][0].vx_b = hp[i][1].vx_b;
	  hp[i][nx].vx_b = hp[i][nx-1].vx_b;
	}
#pragma omp for schedule(static)
	for (i=1;i<ny;i++) {
	  for (j=0;j<nx;j++) {
	    ty = .5*(cells[i-1][j].tby+cells[i][j].tby);
	    ts = .5*(cells[i-1][j].ts+cells[i][j].ts);
	    vb = .5*(cells[i-1][j].vb+cells[i][j].vb);
	    if (vp[i][j].ice > 1.0) vp[i][j].vy_b = (1.0-vbfac)*vp[i][j].vy_b + vbfac*vb*ty/(ts+1.0e-6);
	    else vp[i][j].vy_b = 0.0;
	  }
	  vp[0][i].vy_b = vp[1][i].vy_b;
	  vp[ny][i].vy_b = vp[ny-1][i].vy_b;
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

void smooth_snow_field(int ii,int jj,celltype **cells,double maxslope)
{

  int i,ni,nj;
  double dd,gift,Hd,Hr,sslope;
  int ne = cells[ii][jj].nne;

  /*loop neighboring cells*/
  for (i=0;i<ne;i++) {
    ni = cells[ii][jj].ne_i[i];
    nj = cells[ii][jj].ne_j[i];
    dd = cells[ii][jj].ndist[i];
    Hd = cells[ii][jj].topice+cells[ii][jj].snow;
    Hr = cells[ni][nj].topice+cells[ni][nj].snow;
    sslope = (Hd-Hr)/dd;
    if (sslope > maxslope) {
      gift = .5*(Hd-Hr-dd*maxslope);
      if (gift < 0.0) gift = 0.0;
      if (gift > cells[ii][jj].snow) gift = cells[ii][jj].snow;
      cells[ni][nj].snow += gift;
      cells[ii][jj].snow -= gift;
      if (gift > 1.0e-3) smooth_snow_field(ni,nj,cells,maxslope);
    }/*if*/
  }/*i*/

}

double get_mrate(mproptype mprop,double h)
{

  /*interpolate mass balance function*/
  int i;
  double fac,mrate;

  i = 0;
  while ((h > mprop.Mrate[0][i])&&(i < mprop.nMrate-1)) i += 1;
  if ((i == 0)||(i == mprop.nMrate-1)) mrate = mprop.Mrate[1][i];
  else {
    fac = (h - mprop.Mrate[0][i-1])/(mprop.Mrate[0][i]-mprop.Mrate[0][i-1]);
    mrate = (1.0-fac)*mprop.Mrate[1][i-1] + fac*mprop.Mrate[1][i];
  }

  return(mrate);


}

void ice_accumulation_rate(celltype **cells,meshtype mesh,mproptype mprop)
{

  int i,j;
  
  int nx = mesh.nx;
  int ny = mesh.ny;
  double accrate,topofac,aprop,T;
  double sumaprop = 0.0;
  double sumtprop = 0.0;

  int mtype = mprop.mtype;
  int doava = mesh.doavalance;
  double lrate = mprop.lrate;
  double T0 = mprop.T0;
  double Tsl = mprop.Tsl;
  double maxacc = mprop.maxacc;
  double accgrad = mprop.accgrad;
  double avaslope = mprop.avaslope;
  double avacurv = mprop.avacurv;
  
  /*printf("doava = %d\n",doava);*/


  /*loop cells to determine accumulation rates not corrected for topography*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if (mtype == 1) {
	accrate = get_mrate(mprop,cells[i][j].topice);
	if (accrate < 0.0) accrate = 0.0;
      }
      else if (mtype == 2) {
	T = cells[i][j].Tair;
	if (T < Tsl) accrate = -accgrad*(T-Tsl);
	else accrate = 0.0;
      }
      else {
	accrate = cells[i][j].mrate;
	if (accrate < 0.0) accrate = 0.0;
      }
      if (accrate > cells[i][j].maxacc) accrate = cells[i][j].maxacc;
      cells[i][j].accrate = accrate;


      if (doava > 0) {
	if (accrate > 0.0) {
	  sumaprop += cells[i][j].accrate;
	  if (cells[i][j].lee >= 2.0) sumtprop += cells[i][j].accrate;
	}
      }

    }/*j*/
  }/*i*/

  if (doava > 0) {

    /*accumulation boost due to topography*/
    topofac = (sumaprop+1.0e-6)/(sumtprop+1.0e-6);

    /*loop cells to correct for topography*/
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i][j].lee < 2.0) cells[i][j].accrate = 0.0;
	else {
	  cells[i][j].tprop = topofac;
	  cells[i][j].accrate *= cells[i][j].tprop;
	}
	  

      }
    }

  }

}


void ice_accumulation_rate2(celltype **cells,meshtype mesh,iproptype iprop,mproptype mprop,long **cascade)
{

  int i,j;
  int ii,jj,ri,rj;
  double T,newsnow;

  int nx = mesh.nx;
  int ny = mesh.ny;
  int mtype = mprop.mtype;
  double lrate = mprop.lrate;
  double T0 = mprop.T0;
  double Tsl = mprop.Tsl;
  double maxacc = mprop.maxacc;
  double accgrad = mprop.accgrad;
  double avaslope = mprop.avaslope;
  double avacurv = mprop.avacurv;
  double maxslope = mprop.maxslope;

  double dt = 1.0; /*dt=1 is used for computing acc rates*/

#pragma omp parallel shared(cells,mprop) private(i,j,T,newsnow) firstprivate(T0,lrate,maxacc,accgrad,maxslope,avaslope,mtype,dt)
  {


    double T = 0.0;

    /*compute acuumulation*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].snow = 0.0;
	if (mtype == 1) {
	  newsnow = get_mrate(mprop,cells[i][j].topice);
	  if (newsnow < 0.0) newsnow = 0.0;
	  newsnow *= dt;
	}
	else if (mtype == 2) {
	  T = cells[i][j].Tair;
	  if (T < Tsl) newsnow = -accgrad*(T-Tsl)*dt;
	  else newsnow = 0.0;
	}
	else {
	  newsnow = cells[i][j].mrate;
	  if (newsnow < 0.0) newsnow = 0.0;
	  newsnow *= dt;
	}
	if (newsnow > cells[i][j].maxacc*dt) newsnow = cells[i][j].maxacc*dt;
	cells[i][j].snow += newsnow;


	/*if ((cells[i][j].tslope > avaslope)||((cells[i][j].curv < avacurv)&&(cells[i][j].ice < 10.0))) cells[i][j].snow = 0.0;*/
	/*if ((cells[i][j].streamorder < 5)&&(cells[i][j].ice < 10.0)) cells[i][j].snow = 0.0;*/



      }/*j*/
    }/*i*/


  }/*pragma*/



  /*do avalances*/
  if (mesh.doavalance > 0) {

    /*avalance snow downhill*/
    for (i=0;i<nx*ny;i++) {
      ii = cascade[i][0];
      jj = cascade[i][1];
      /*if ((cells[ii][jj].tslope > avaslope)||((cells[ii][jj].streamorder < 5)&&(cells[ii][jj].ice < 10.0))) {*/
      if ((cells[ii][jj].tslope > avaslope)||(cells[ii][jj].streamorder < 5)) {
	ri = cells[ii][jj].ri;
	rj = cells[ii][jj].rj;
	if (cells[ri][rj].snow > 0.0) cells[ri][rj].snow += cells[ii][jj].snow; /*only avalance to accumulation zones*/
	cells[ii][jj].snow = 0.0;
      }
    }
    
    /*snooth the snowfield to critical slope*/
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i][j].snow > 1e-3) smooth_snow_field(i,j,cells,maxslope);
      }
    }

  }

  /*save the accumulation rate*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].accrate = cells[i][j].snow/dt;
      if (cells[i][j].accrate > 2.0*cells[i][j].maxacc) cells[i][j].accrate = 2.0*cells[i][j].maxacc;
    }
  }


}

void ice_accumulation_rate3(celltype **cells,meshtype mesh,mproptype mprop)
{

  int i,j;
  
  int nx = mesh.nx;
  int ny = mesh.ny;
  double accrate,topofac,aprop,T;
  double sumaprop = 0.0;
  double sumtprop = 0.0;

  int mtype = mprop.mtype;
  int doava = mesh.doavalance;
  double lrate = mprop.lrate;
  double T0 = mprop.T0;
  double Tsl = mprop.Tsl;
  double maxacc = mprop.maxacc;
  double accgrad = mprop.accgrad;
  double avaslope = mprop.avaslope;
  double avacurv = mprop.avacurv;
  
  /*printf("doava = %d\n",doava);*/


  /*loop cells to determine accumulation rates not corrected for topography*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if (mtype == 1) {
	accrate = get_mrate(mprop,cells[i][j].topice);
	if (accrate < 0.0) accrate = 0.0;
      }
      else if (mtype == 2) {
	T = cells[i][j].Tair;
	if (T < Tsl) accrate = -accgrad*(T-Tsl);
	else accrate = 0.0;
      }
      else {
	accrate = cells[i][j].mrate;
	if (accrate < 0.0) accrate = 0.0;
      }
      if (accrate > cells[i][j].maxacc) accrate = cells[i][j].maxacc;
      cells[i][j].accrate = accrate;


      if (doava > 0) {
	if (accrate > 0.0) {
	  sumaprop += cells[i][j].accrate;
	  if (cells[i][j].tslope < avaslope) sumtprop += cells[i][j].accrate;
	}
      }

    }/*j*/
  }/*i*/

  if (doava > 0) {

    /*accumulation boost due to topography*/
    topofac = (sumaprop+1.0e-6)/(sumtprop+1.0e-6);

    /*loop cells to correct for topography*/
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i][j].tslope >= avaslope) cells[i][j].accrate = 0.0;
	else {
	  cells[i][j].tprop = topofac;
	  /*cells[i][j].accrate *= cells[i][j].tprop;*/
	}
	  

      }
    }

  }

}


void mass_balance(celltype **cells,meshtype *mesh,iproptype iprop,mproptype mprop,hwproptype hwprop,double dt)
{

  int i,j;
  double T,melt,accumulation,smelt,bmelt;
  double meanice_b,dice = 0.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int doglacialsedi = (*mesh).doglacialsedi;
  int dodebrisablation = (*mesh).dodebrisablation;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int mtype = mprop.mtype;
  double g = (*mesh).gravity;
  double lrate = mprop.lrate;
  double T0 = mprop.T0;
  double Tsl = mprop.Tsl;
  double ablgrad = mprop.ablgrad; 
  double maxabla = mprop.maxabla;
  double L = iprop.latentheat; 
  double sedifrac = mprop.sedifrac; 
  double Ldebris = mprop.Ldebris;

#pragma omp parallel shared(cells,mesh,mprop,hwprop,dice,meanice_b) private(i,j,T,melt,smelt,bmelt,accumulation) firstprivate(nx,ny,lrate,T0,ablgrad,maxabla,dt,g,L,mtype,dx,dy,doglacialsedi,sedifrac,dodebrisablation,Ldebris)
  {

    /*compute acuumulation and ablation*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	
	/*ice accumulation*/
	accumulation = dt*cells[i][j].accrate;

	if (mtype == 1) {

	  smelt = get_mrate(mprop,cells[i][j].topice);
	  if (smelt > 0.0) smelt = 0.0;
	  smelt *= -1.0;
	  if (smelt > maxabla) smelt = maxabla;

	}
	else if (mtype == 2) {

	  /*Atmosphere temperature*/
	  T = cells[i][j].Tair;

	  /*ablation*/
	  if (cells[i][j].maxacc == 0.0) smelt = maxabla;
	  else {
	    if (T > Tsl) smelt = ablgrad*(T-Tsl);
	    else smelt = 0.0;
	    if (smelt > maxabla) smelt = maxabla;
	  }
	}
	else {

	  smelt = cells[i][j].mrate;
	  if (smelt > 0.0) smelt = 0.0;
	  smelt *= -1.0;
	  if (smelt > maxabla) smelt = maxabla;

	}

	/*basal melting*/
	bmelt = -cells[i][j].Mb; /*Mb is negative in thermal.c*/

	/*reduce surface melt due to debris cover*/
	if (dodebrisablation) smelt *= Ldebris/(cells[i][j].Vs[0]+Ldebris);

	/*scale melt*/
	smelt *= dt;
	bmelt *= dt;

	/*total melt*/
	melt = smelt + bmelt;

	/*limit melting*/
	if ((melt > (cells[i][j].ice+accumulation))&&(melt > 0.0)) {
	  bmelt *= (cells[i][j].ice+accumulation)/melt;
	  smelt *= (cells[i][j].ice+accumulation)/melt;
	  melt = cells[i][j].ice + accumulation;
	}

	/*store melt rate*/
	cells[i][j].meltrate = melt/dt;
	cells[i][j].Ms = (-smelt + accumulation)/dt;

	/*store input to glacial hydrology*/
	/*if ((*mesh).doglacialhydrology > 0) cells[i][j].dHw = hwprop.a2w*smelt/dt + bmelt/dt;
	  else cells[i][j].dHw = 0.0;*/

	/*do mass balance*/
	cells[i][j].ice += accumulation - melt;
	if (doglacialsedi > 0) cells[i][j].Vs[0] += sedifrac*accumulation;

	/*pass melts to hydrology*/
	/*if ((*mesh).doglacialhydrology > 0) cells[i][j].dHw = (hwprop.a2w*smelt + bmelt)/dt;*/
	if ((*mesh).doglacialhydrology > 0) cells[i][j].Hgw += hwprop.tscale*(hwprop.a2w*smelt + bmelt);

      }/*j*/
    }/*i*/


    /*compute total ice volume change*/
#pragma omp for schedule(static) reduction(+:dice)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*dice += cells[i][j].meltrate;*/
      dice += cells[i][j].accrate-cells[i][j].meltrate;
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
  double kv = 0.0;

#pragma omp parallel shared(cells,vp,hp,mesh) private(i,j,dH) firstprivate(nx,ny,dx,dy,dt,kv)
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
    for (j=1;j<nx;j++) {
      dH = hp[i][j].vx*hp[i][j].ice*dt/dx;
      if (cells[i][j].topice < cells[i][j-1].topice) dH -= kv*(cells[i][j].ice-cells[i][j-1].ice);
      /*hdiff = .5*(cells[i][j-1].topice-cells[i][j].topice);
      if ((dH > 0.0)&&(dH > hdiff)) { dH = hdiff; if (dH < 0.0) dH = 0.0;}
	if ((dH < 0.0)&&(dH < hdiff)) { dH = hdiff; if (dH > 0.0) dH = 0.0;}
      if ((dH > 0.0)&&(cells[i][j-1].ice < 0.0)) dH = 0.0;
      if ((dH < 0.0)&&(cells[i][j].ice < 0.0)) dH = 0.0;*/
      if (dH < 0.0) {
	if (cells[i][j].ice <= 0.0) dH = 0.0;
	else if (dH < -cells[i][j].ice) dH = -cells[i][j].ice;
      }
      else {
	if (cells[i][j-1].ice <= 0.0) dH = 0.0;
	else if (dH > cells[i][j-1].ice) dH = cells[i][j-1].ice;
	}
      hp[i][j].dH = dH;
    }
  }

  /*loop v-points*/
#pragma omp for schedule(static)
  for (i=1;i<ny;i++) {
    for (j=0;j<nx;j++) {
      dH = vp[i][j].vy*vp[i][j].ice*dt/dy;
      if (cells[i][j].topice < cells[i-1][j].topice) dH -= kv*(cells[i][j].ice-cells[i-1][j].ice);
      /*hdiff = .5*(cells[i-1][j].topice-cells[i][j].topice);
      if ((dH > 0.0)&&(dH > hdiff)) { dH = hdiff; if (dH < 0.0) dH = 0.0;}
      if ((dH < 0.0)&&(dH < hdiff)) { dH = hdiff; if (dH > 0.0) dH = 0.0;}
      if ((dH > 0.0)&&(cells[i-1][j].ice < 0.0)) dH = 0.0;
      if ((dH < 0.0)&&(cells[i][j].ice < 0.0)) dH = 0.0;*/
      if (dH < 0.0) {
	if (cells[i][j].ice <= 0.0) dH = 0.0; 
	else if (dH < -cells[i][j].ice) dH = -cells[i][j].ice;
      }
      else {
	if (cells[i-1][j].ice <= 0.0) dH = 0.0;
	else if (dH > cells[i-1][j].ice) dH = cells[i-1][j].ice;
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

void glacial_erosion(celltype **cells,hptype **hp,vptype **vp,iproptype iprop,meshtype *mesh,double dt)
{

  int i,j;
  double ub,slope,pe,qtot=0.0,atot=0.0,meanqrate=0.0,meanarate=0.0;
  double mf = iprop.mf;
  double qfac = iprop.qfac;
  double Ka = iprop.Ka;
  double ap = iprop.ap;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
 

#pragma omp parallel shared(cells,mesh,meanqrate,meanarate) private(i,j,ub,slope,pe) firstprivate(nx,ny,mf,Ka,dt,qtot,atot)
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
	if ((cells[i][j].sliding > 1.0)&&(cells[i][j].margin <= 0)) {
	  ub = cells[i][j].sliding;
	  slope = -(cells[i][j].dbdx*cells[i][j].vx_b+cells[i][j].dbdy*cells[i][j].vy_b)/ub;
	  pe = 0.01*cells[i][j].tn;
	  if (pe < 0.0) pe = 0.0; if (pe > 10.0) pe = 10.0;
	  if (slope > -0.55) {
	    if (slope > 1.0) slope = 1.0;
	    if (mf == 1.5) {
	      qtot = 2.8e-2*pow(pe,2.0)*pow(ub,0.7)*pow(slope+0.55,1.0);
	      /*if (pe < 0.7) qtot = 0.11*pow(pe,3.5)*pow(ub,0.35)*pow(slope+0.55,3.9);
		else qtot = 5.4e-3*pow(pe,1.4)*pow(ub,0.97)*pow(slope+0.55,3.9);*/
	    }
	    else if (mf == 2.0) {
	      if (pe < 0.7) qtot = 0.018*pow(pe,3.7)*pow(ub,0.42)*pow(slope+0.55,3.8);
	      else qtot = 0.001*pow(pe,1.9)*pow(ub,1.0)*pow(slope+0.55,4.0);
	    }
	    else if (mf == 2.25) {
	      if (pe < 0.7) qtot = 0.003*pow(pe,3.93)*pow(ub,0.5)*pow(slope+0.55,4.2);
	      else qtot = 0.0002*pow(pe,2.4)*pow(ub,1.0)*pow(slope+0.55,3.8);
	    }
	    else if (mf == 2.5) {
	      if (pe < 0.7) qtot = 0.003*pow(pe,4.0)*pow(ub,0.5)*pow(slope+0.55,3.8);
	      else qtot = 0.0002*pow(pe,2.4)*pow(ub,1.0)*pow(slope+0.55,3.9);
	    }
	    else if (mf == 3.0) {
	      /*if (pe < 0.7) qtot = 3.6e-4*pow(pe,3.9)*pow(ub,0.66)*pow(slope+0.55,4.4);
		else qtot = 3.3e-5*pow(pe,2.9)*pow(ub,1.0)*pow(slope+0.55,3.8);*/
	      /*qtot = 3.7e-4*pow(pe,3.0)*pow(ub,0.9)*pow(slope+0.55,2.1);*/
	      qtot = 1.0e-4*pow(pe,3.0)*pow(ub,1.0)*pow(slope+0.55,2.0);
	    }
	    else if (mf == 5.0) {
	      /*if (pe < 0.7) qtot = 9.5e-7*pow(pe,3.8)*pow(ub,1.0)*pow(slope+0.55,5.5);
		else qtot = 6.0e-8*pow(pe,4.3)*pow(ub,1.2)*pow(slope+0.55,4.0); */
	      qtot = 1.0e-6*pow(pe,3.6)*pow(ub,1.2)*pow(slope+0.55,2.6);
	    }
	    else if (mf == 10.0) {
	      qtot = 5.0e-6*pow(pe,4.0)*pow(ub,0.7)*pow(slope+0.55,2.0);
	    }
	    else {
	      /*printf("In glacial erosion: value of mf not supported\n");*/
	      qtot = 0.0;
	    }
	    qtot *= qfac*sqrt(1.0+slope*slope);
	  }
	  else qtot = 0.0;

	  /*printf("ub = %4.4e, pe = %4.4e, slope = %4.4e, qtot = %4.4e\n",ub,pe,slope,qtot);*/

	  /*abrasion*/
	  atot = Ka*pow(ub,ap)*sqrt(1.0+slope*slope);
	  /*if (slope > -0.55) qtot = 2.4e-6*pow(pe,3.2)*pow(ub,0.8)*pow(slope+0.55,4.0); else qtot = 0.0;*/

	  /*remove sediment*/
	  /*
	  cells[i][j].sedi -= atot*dt;
	  if (cells[i][j].sedi < 0.0) {
	    atot = -cells[i][j].sedi/dt;
	    cells[i][j].sedi = 0.0;
	    }*/

	  /*damp quarrying by sediment cover*/
	  /*qtot *= exp(-cells[i][j].sedi/0.5);
	  atot *= exp(-cells[i][j].sedi/0.5);
	  */
	  /*erode cell*/
	  cells[i][j].quarrying_rate = qtot;
	  cells[i][j].abrasion_rate = atot;
	  cells[i][j].quarrying += qtot*dt;
	  cells[i][j].abrasion += atot*dt;
	  cells[i][j].bed -= (qtot+atot)*dt;

	  /*cells[i][j].sedi += (qtot+atot)*dt;*/

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

void get_margin(celltype **cells,meshtype mesh)
{

  int i,j;
  double minice = 10.0;
  int nx = mesh.nx;
  int ny = mesh.ny;

  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {

      if (cells[i][j].ice < minice) cells[i][j].margin = 0;
      else {

	if (cells[i][j-1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i+1][j].ice < minice) cells[i][j].margin = 1;
	else if (cells[i][j+1].ice < minice) cells[i][j].margin = 1;
	else if (cells[i-1][j].ice < minice) cells[i][j].margin = 1;
	else cells[i][j].margin = -1;

      }/*else*/

    }/*j*/
  }/*i*/

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


void glacial_sediment_transport(celltype **cells,hptype **hp,vptype **vp,meshtype mesh,double dt)
{

  int i,j,k,t;
  double dz,zz,ice;
  double ms,mb,dc,cs;
  double vx,vy,vz;

  int ny = mesh.ny;
  int nx = mesh.nx;
  double dx = mesh.dx;
  double dy = mesh.dy;
  int dosediment = mesh.dosediment;

  int Nz = 20;
  double minice = 20.0;

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
	      cells[i][j].Vs[k] += dc;
	      cells[i][j].Vs[k-1] -= dc;
	    }
	    else if (vz > 0.0) {
	      dc = dtz*vz*cells[i][j].Vs[k-1]/dz; /*dc>0*/
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
	    cells[i][j].Vs[k] += dc;
	    cells[i][j-1].Vs[k] -= dc;
	  }/*k*/
	}/*if*/
	else if (vx > 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vx*cells[i][j-1].Vs[k]/dx;/*dc>0*/
	    cells[i][j-1].Vs[k] -= dc;
	    cells[i][j].Vs[k] += dc;
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
	    cells[i][j].Vs[k] += dc;
	    cells[i-1][j].Vs[k] -= dc;
	  }/*k*/
	}/*if*/
	else if (vy > 0.0) {
	  for (k=0;k<Nz;k++) {
	    dc = dt*vy*cells[i-1][j].Vs[k]/dy; /*dc>0*/
	    cells[i-1][j].Vs[k] -= dc;
	    cells[i][j].Vs[k] += dc;
	  }/*k*/
	}/*if*/
      }/*if*/
    }/*j*/
  }/*i*/
  
  /*transfer sediment to sedi array*/
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
