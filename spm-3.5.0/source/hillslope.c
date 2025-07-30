


void hillslope_production(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hproptype hprop,double dt)
{

  int i,j;
  double sdiff,fac,ero,maxero;
  double meanerate = 0.0;

  double sc;
  double Ke = hprop.Ke;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int dosediment = (*mesh).dosediment;

#pragma omp parallel shared(cells,hp,vp,meanerate) private(i,j,sdiff,fac,maxero,ero,sc) firstprivate(nx,ny,dx,dy,dt,dosediment,Ke)
  {

    /***************** Erosion (critical slope) *****************/
      /*loop h-points for bedrock erosion*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) { 
	  if (cells[i][j-1].include+cells[i][j].include > 1.0) {

	    sc = 0.5*(cells[i][j-1].phi+cells[i][j].phi);

	    /*if (fabs(hp[i][j].dbdx) > sc) {*/

	      fac = 1.0 - pow(hp[i][j].dbdx/sc,2.0); if (fac < 1.0e-4) fac = 1.0e-4;
	      sdiff = Ke/fac;
	      ero = sdiff*hp[i][j].dbdx*dt/dx*sqrt(1.0+pow(cells[i][j].bslope,2.0));
#pragma omp critical
	      {
		if (ero < 0.0) {
		  /*maxero = cells[i][j-1].bed-cells[i][j].bed-sc*dx;
		  if (maxero < 0.0) maxero = 0.0;
		  if (ero < -maxero) ero = -maxero;*/
		  cells[i][j-1].bed += ero;
		  cells[i][j-1].hillslope_erosion -= ero;
		  cells[i][j-1].hillslope_erate = -ero/dt;
		  if (dosediment > 0) cells[i][j-1].sedi -= ero;
		}
		else {
		  maxero = cells[i][j].bed-cells[i][j-1].bed-sc*dx;
		  if (maxero < 0.0) maxero = 0.0;
		  if (ero > maxero) ero = maxero;
		  cells[i][j].bed -= ero;
		  cells[i][j].hillslope_erosion += ero;
		  cells[i][j].hillslope_erate = ero/dt;
		  if (dosediment > 0) cells[i][j].sedi += ero;
		}
	      } 
	      /*}*/
	  }
	}
      }


    /*loop v-points*/
#pragma omp for schedule(static)
      for (j=0;j<nx;j++) {
	for (i=1;i<ny;i++) {
	  if (cells[i-1][j].include+cells[i][j].include > 1.0) {

	    sc = 0.5*(cells[i-1][j].phi+cells[i][j].phi);

	    /*if (fabs(vp[i][j].dbdy) > sc) {*/


	      fac = 1.0 - pow(vp[i][j].dbdy/sc,2.0); if (fac < 1.0e-4) fac = 1.0e-4;
	      sdiff = Ke/fac;
	      ero = sdiff*vp[i][j].dbdy*dt/dy*sqrt(1.0+pow(cells[i][j].bslope,2.0));
#pragma omp critical
	      {
		if (ero < 0.0) {
		  /*maxero = cells[i-1][j].bed-cells[i][j].bed-sc*dy;
		  if (maxero < 0.0) maxero = 0.0;
		  if (ero < -maxero) ero = -maxero;*/
		  cells[i-1][j].bed += ero;
		  cells[i-1][j].hillslope_erosion -= ero;
		  cells[i-1][j].hillslope_erate = -ero/dt;
		  if (dosediment > 0) cells[i-1][j].sedi -= ero;
		}
		else {
		  maxero = cells[i][j].bed-cells[i-1][j].bed-sc*dy;
		  if (maxero < 0.0) maxero = 0.0;
		  if (ero > maxero) ero = maxero;
		  cells[i][j].bed -= ero;
		  cells[i][j].hillslope_erosion += ero;
		  cells[i][j].hillslope_erate = ero/dt;
		  if (dosediment > 0) cells[i][j].sedi += ero;
		}
	      } 
	      /*}*/
	  }
	}
      }
    

  /*compute mean erosion rate and transport change*/
#pragma omp for schedule(static) reduction(+:meanerate)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanerate += cells[i][j].hillslope_erate;
    }
  }
#pragma omp single
  {
    meanerate /= (double)((*mesh).nc);
    (*mesh).mean_hillslope_erate = meanerate;

  }

  }

}

void hillslope_transport(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hproptype hprop,double dt)
{

  int i,j;
  double sdiff,dHs,mean_dHs=0.0,fac;

  double Ks = hprop.Ks;
  double sc = hprop.sc;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double minice = 10.0;


  /*********************** sediment transport *********************/


#pragma omp parallel shared(cells,hp,vp,mean_dHs) private(i,j,sdiff,dHs,fac) firstprivate(nx,ny,dx,dy,dt,Ks,sc,minice)
  {


    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dHs = 0.0;
      }
    }
    
    /*loop h-points for horizontal transport*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) {
	  fac = 1.0 - pow(hp[i][j].dtdx/sc,2.0); if (fac < 0.001) fac = 0.001;
	  sdiff = Ks/fac;
	  /*dHs = -.5*sdiff*hp[i][j].dbdx*(cells[i][j-1].sedi+cells[i][j].sedi)*dt/dx;*/
	  dHs = -sdiff*hp[i][j].dtdx*dt/dx;
	  if (dHs < 0.0) {
	    if (cells[i][j].sedi <= 0.0) dHs = 0.0;
	    else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	  }
	  else {
	    if (cells[i][j-1].sedi <= 0.0) dHs = 0.0;
	    else if (dHs > cells[i][j-1].sedi) dHs = cells[i][j-1].sedi;
	  }
#pragma omp critical
	  {
	    cells[i][j].dHs += dHs;
	    cells[i][j-1].dHs -= dHs;
	  }
	}
      }

      /*loop v-points*/
#pragma omp for schedule(static)
      for (j=0;j<nx;j++) {
	for (i=1;i<ny;i++) {
	  fac = 1.0 - pow(vp[i][j].dtdy/sc,2.0); if (fac < 0.001) fac = 0.001;
	  sdiff = Ks/fac;
	  /*dHs = -.5*sdiff*vp[i][j].dbdy*(cells[i-1][j].sedi+cells[i][j].sedi)*dt/dy;*/
	  dHs = -sdiff*vp[i][j].dtdy*dt/dy;
	  if (dHs < 0.0) {
	    if (cells[i][j].sedi <= 0.0) dHs = 0.0; 
	    else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	  }
	  else {
	    if (cells[i-1][j].sedi <= 0.0) dHs = 0.0;
	    else if (dHs > cells[i-1][j].sedi) dHs = cells[i-1][j].sedi;
	  }
#pragma omp critical
	  {
	    cells[i][j].dHs += dHs;
	    cells[i-1][j].dHs -= dHs;
	  }
	}
      }
      
    
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  cells[i][j].sedi += cells[i][j].dHs;
	}
      }
    
  /*compute mean erosion rate and transport change*/
#pragma omp for schedule(static) reduction(+:mean_dHs) 
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      mean_dHs += fabs(cells[i][j].dHs);
    }
  }
#pragma omp single
  {
    mean_dHs /= (double)((*mesh).nc)*dt*2.0;
    (*mesh).mean_dHs_hillslope = mean_dHs;
  }


  }

}


void nslide(int ti,int tj,celltype **cells,long *nsc,long **slist,double x0,double y0,double h0,double *maxH,double *maxbeta)
{

  long i;
  int ni,nj;
  double hp,dist,Hc,beta; 
  double xc = cells[ti][tj].x;
  double yc = cells[ti][tj].y;
  double hc = cells[ti][tj].topsedi;

  /*printf("ti = %d, tj = %d, nne = %d\n",ti,tj,cells[ti][tj].nne);*/

  /*loop neighbours*/
  for (i=0;i<cells[ti][tj].nne;i++) {

    /*identify neighbour*/
    ni = cells[ti][tj].ne_i[i];
    nj = cells[ti][tj].ne_j[i];

    /*printf("ni = %d, nj = %d\n",ni,nj);*/


    /*if not allready in slide*/
    if (cells[ni][nj].inslide == 0) {
    
      /*distance to neigbour*/
      dist = sqrt(pow(cells[ni][nj].x-xc,2.0)+pow(cells[ni][nj].y-yc,2.0));

      /*Critical hillslope height*/
      Hc = hc + (double)tan(cells[ni][nj].phi*pi/180.0)*dist;
      /*Hc = hc + (double)tan(phi0)*dist;*/

      /*printf("Hc = %4.4e, h0 = %4.4e ",Hc,h0);*/

      /*Height of cell*/
      hp = cells[ni][nj].topsedi;


      /*printf("hp = %4.4e, Hc = %4.4e\n",hp,Hc);*/

      /*if potentially unstable*/
      if (hp > Hc) {

	/*registre cell*/
	slist[(*nsc)][0] = ni;
	slist[(*nsc)][1] = nj;
	cells[ni][nj].inslide = 1;
	cells[ni][nj].rfail = hp-Hc;

	/*add to number*/
	(*nsc) += 1;
	/*printf("nsc = %ld\n",(*nsc));*/

	/*increase hillslope height*/
	if ((*maxH) < hp) (*maxH) = hp;

	/*distance to toe*/
	dist = sqrt(pow(cells[ni][nj].x-x0,2.0)+pow(cells[ni][nj].y-y0,2.0));

	/*increase hillslope angle*/
	beta = (double)atan((hp-h0)/dist);
	if ((*maxbeta) < beta) (*maxbeta) = beta;

	/*recursive call*/
	nslide(ni,nj,cells,nsc,slist,x0,y0,h0,maxH,maxbeta);
	
      }/*if*/

      /*printf("\n");*/

    }/*if*/

  }/*i*/

}




void landslide(celltype **cells,meshtype *mesh,hproptype hprop,long **slist,double dt)
{

  char file[200];
  long i,j,k;
  int ti,tj,ni,nj;
  long nsc,nf;
  double HH,Hc,dH,dHmax,maxbeta,theta,phi,Pf,lfrac,ran,mdist;
  double maxH,x0,y0,h0;
  double Atot,Vtot,jgift,newsedi,Vmean;
  double ero;

  double gamma = hprop.gamma;
  int Nc = hprop.Nc;
  double kt = 1.0e-5;
  FILE *fg;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double area = dx*dy;

  int dosediment = 0;

  /*initialize*/
  for (i=0;i<ny;i++) 
    for (j=0;j<nx;j++) {
      cells[i][j].inslide = 0;
      cells[i][j].rfail = 0.0;
    }

  /*mean landslide volume*/
  Vmean = 0.0;

  /*loop some random cells*/
  for (i=1;i<=Nc;i++) {

    /*random double between 0 and 1*/
    ran = (double)rand()/(double)RAND_MAX;

    /*random integer between 1 and ny*/
    ti = (long)floor(ran*(double)ny);

    /*random double between 0 and 1*/
    ran = (double)rand()/(double)RAND_MAX;

    /*random integer between 1 and nx*/
    tj = (long)floor(ran*(double)nx);

    /*if not allready in slide*/
    if ((ti >=0)&&(ti<ny)&&(tj>=0)&&(tj<nx)&&cells[ti][tj].inslide == 0) {

      nsc = 0;
      x0 = cells[ti][tj].x;
      y0 = cells[ti][tj].y;
      h0 = cells[ti][tj].topsedi;
      maxH = h0;
      maxbeta = 0.0;
      phi = cells[ti][tj].phi*pi/180.0;

      /*recursive call*/
      nslide(ti,tj,cells,&nsc,slist,x0,y0,h0,&maxH,&maxbeta);

      /*printf("nsc = %ld\n",nsc);*/


      /*if potential for slide*/
      if (nsc > 0) {
	
	/*Hillslope height*/
	HH = maxH - h0;

	/*Critical height*/
	Hc = gamma*sin(maxbeta)*cos(phi)/(1.0-cos(maxbeta-phi));

	/*Propability*/
	Pf = HH/Hc;
	/*Pf = HH/Hc + kt*(time - lastslide[tp]);*/

	/*pick random number between 0 and 1*/
	ran = (double)rand()/(double)RAND_MAX;

	/*if failure*/
	if (ran < Pf) {
	  
	  /*critical angle*/
	  /*theta = 0.5*(maxbeta + phi0);*/

	  Atot = 0.0;
	  Vtot = 0.0;

	  /*loop unstable cells and collect material*/
	  for (j=0;j<nsc;j++) {
	    
	    /*identify cell*/
	    ni = slist[j][0];
	    nj = slist[j][1];

	    /*distance to cell*/
	    /*mdist = sqrt(pow(cells[ni][nj].x-x0,2.0)+pow(cells[ni][nj].y-y0,2.0));*/

	    /*Excess height*/
	    /*dH = cells[ni][nj].topsedi - (h0+(double)tan(theta)*mdist);*/
	    dH = cells[ni][nj].rfail;

	    /*erosion*/
	    if (dH > 0.1) {
	      
	      if (dH > cells[ni][nj].sedi) {
		newsedi = dH - cells[ni][nj].sedi;
		cells[ni][nj].bed -= newsedi;
		cells[ni][nj].landslide_erosion += newsedi;
		if (dosediment > 0) {
		  cells[ni][nj].sedi += newsedi;
		}
	      }

	      Atot += area;
	      Vtot += area*dH;
	      /*printf("nsc = %ld,tn = %ld, dH = %4.4e\n",nsc,tn,dH);*/

	    }/*if erosion*/

	  }/*j*/

	  /*report slide to file*/
	  /*
	  sprintf(file,"%s/output/landslides.dat",path);
	  if ((fg = fopen(file,"at")) == NULL) {
	    nowrite = 1;
	  }
	  else {
	    fprintf(fg,"%2.5e %2.5e %2.5e %ld %ld\n",time,Atot,Vtot,tp,nsc);
	    fclose(fg);
	    }
	  */
	  /*report big landslide*/
	  /*
	  if (Vtot > 5.0e8) {

	    sprintf(file,"%s/output/bigslides.dat",path);
	    if ((fg = fopen(file,"at")) == NULL) {
	      nowrite = 1;
	    }
	    else {
	      fprintf(fg,"%2.5e %2.5e %2.5e %ld ",time,Atot,Vtot,nsc);
	      for (j=1;j<=nsc;j++) fprintf(fg,"%d ",slist[j]);
	      fprintf(fg,"\n");
	      fclose(fg);
	    }

	  }*/

	  Vmean += Vtot;


      }/*if failure*/

      }/*if nsc > 0*/

    }/*if not in slide*/

  }/*i*/

  (*mesh).mean_landslide_erate = Vmean/((*mesh).L*(*mesh).H*dt);

}/*landslide*/


void weathering(celltype **cells,meshtype *mesh,hproptype hprop,double dt)
{

  int i,j;

  double Kw = hprop.Kw;
  double Ls = hprop.Ls;
  double ero,meanwrate = 0.0;
  double minice = 10.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int dosediment = (*mesh).dosediment;


#pragma omp parallel shared(cells,meanwrate) private(i,j,ero) firstprivate(nx,ny,dt,Kw,Ls,minice,dosediment)
  {

#pragma omp for schedule(static) reduction(+:meanwrate)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  if ((cells[i][j].ice < minice)&&(cells[i][j].include > 0)) {
	    ero = Kw*exp(-cells[i][j].sedi/Ls)*dt*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	    cells[i][j].bed -= ero;
	    cells[i][j].weathering += ero;
	    cells[i][j].weathering_rate = ero/dt;
	    if (dosediment > 0) cells[i][j].sedi += ero;
	    meanwrate += ero/dt;
	  }
	}
      }
 



#pragma omp single
  {
    meanwrate /= (double)((*mesh).nc);
    (*mesh).mean_weatheringrate = meanwrate;

  }

  } 

}


