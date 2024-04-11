
void initiate_drainage_info(celltype **cells,meshtype mesh)
{

  int i,j;
  int ny = mesh.ny;
  int nx = mesh.nx;
  double dy = mesh.dy;
  double dx = mesh.dx;
  double dd = sqrt(dx*dx+dy*dy);

  /*corners*/

  /*lower left*/
  i = 0; j = 0;
  cells[i][j].nne = 3;
  cells[i][j].ne_i[0] = i+1;
  cells[i][j].ne_i[1] = i+1;
  cells[i][j].ne_i[2] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j+1;
  cells[i][j].ne_j[2] = j+1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dd;
  cells[i][j].ndist[2] = dx;

  /*upper left*/
  i = ny-1; j = 0;
  cells[i][j].nne = 3;
  cells[i][j].ne_i[0] = i-1;
  cells[i][j].ne_i[1] = i-1;
  cells[i][j].ne_i[2] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j+1;
  cells[i][j].ne_j[2] = j+1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dd;
  cells[i][j].ndist[2] = dx;

  /*upper right*/
  i = ny-1; j = nx-1;
  cells[i][j].nne = 3;
  cells[i][j].ne_i[0] = i-1;
  cells[i][j].ne_i[1] = i-1;
  cells[i][j].ne_i[2] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j-1;
  cells[i][j].ne_j[2] = j-1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dd;
  cells[i][j].ndist[2] = dx;

  /*lower right*/
  i = 0; j = nx-1;
  cells[i][j].nne = 3;
  cells[i][j].ne_i[0] = i+1;
  cells[i][j].ne_i[1] = i+1;
  cells[i][j].ne_i[2] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j-1;
  cells[i][j].ne_j[2] = j-1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dd;
  cells[i][j].ndist[2] = dx;

  /*left side*/
  j = 0;
  for (i=1;i<ny-1;i++) {
    cells[i][j].nne = 5;
    cells[i][j].ne_i[0] = i-1;
    cells[i][j].ne_i[1] = i-1;
    cells[i][j].ne_i[2] = i;
    cells[i][j].ne_i[3] = i+1;
    cells[i][j].ne_i[4] = i+1;
    cells[i][j].ne_j[0] = j;
    cells[i][j].ne_j[1] = j+1;
    cells[i][j].ne_j[2] = j+1;
    cells[i][j].ne_j[3] = j+1;
    cells[i][j].ne_j[4] = j;
    cells[i][j].ndist[0] = dy;
    cells[i][j].ndist[1] = dd;
    cells[i][j].ndist[2] = dx;
    cells[i][j].ndist[3] = dd;
    cells[i][j].ndist[4] = dy;
  }

  /*right side*/
  j = nx-1;
  for (i=1;i<ny-1;i++) {
    cells[i][j].nne = 5;
    cells[i][j].ne_i[0] = i-1;
    cells[i][j].ne_i[1] = i-1;
    cells[i][j].ne_i[2] = i;
    cells[i][j].ne_i[3] = i+1;
    cells[i][j].ne_i[4] = i+1;
    cells[i][j].ne_j[0] = j;
    cells[i][j].ne_j[1] = j-1;
    cells[i][j].ne_j[2] = j-1;
    cells[i][j].ne_j[3] = j-1;
    cells[i][j].ne_j[4] = j;
    cells[i][j].ndist[0] = dy;
    cells[i][j].ndist[1] = dd;
    cells[i][j].ndist[2] = dx;
    cells[i][j].ndist[3] = dd;
    cells[i][j].ndist[4] = dy;
  }

  /*lower side*/
  i = 0;
  for (j=1;j<nx-1;j++) {
    cells[i][j].nne = 5;
    cells[i][j].ne_i[0] = i;
    cells[i][j].ne_i[1] = i+1;
    cells[i][j].ne_i[2] = i+1;
    cells[i][j].ne_i[3] = i+1;
    cells[i][j].ne_i[4] = i;
    cells[i][j].ne_j[0] = j-1;
    cells[i][j].ne_j[1] = j-1;
    cells[i][j].ne_j[2] = j;
    cells[i][j].ne_j[3] = j+1;
    cells[i][j].ne_j[4] = j+1;
    cells[i][j].ndist[0] = dx;
    cells[i][j].ndist[1] = dd;
    cells[i][j].ndist[2] = dy;
    cells[i][j].ndist[3] = dd;
    cells[i][j].ndist[4] = dx;
  }

  /*upper side*/
  i = ny-1;
  for (j=1;j<nx-1;j++) {
    cells[i][j].nne = 5;
    cells[i][j].ne_i[0] = i;
    cells[i][j].ne_i[1] = i-1;
    cells[i][j].ne_i[2] = i-1;
    cells[i][j].ne_i[3] = i-1;
    cells[i][j].ne_i[4] = i;
    cells[i][j].ne_j[0] = j-1;
    cells[i][j].ne_j[1] = j-1;
    cells[i][j].ne_j[2] = j;
    cells[i][j].ne_j[3] = j+1;
    cells[i][j].ne_j[4] = j+1;
    cells[i][j].ndist[0] = dx;
    cells[i][j].ndist[1] = dd;
    cells[i][j].ndist[2] = dy;
    cells[i][j].ndist[3] = dd;
    cells[i][j].ndist[4] = dx;
  }

  /*middle*/
  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {
      cells[i][j].nne = 8;
      cells[i][j].ne_i[0] = i;
      cells[i][j].ne_i[1] = i+1;
      cells[i][j].ne_i[2] = i+1;
      cells[i][j].ne_i[3] = i+1;
      cells[i][j].ne_i[4] = i;
      cells[i][j].ne_i[5] = i-1;
      cells[i][j].ne_i[6] = i-1;
      cells[i][j].ne_i[7] = i-1;
      cells[i][j].ne_j[0] = j-1;
      cells[i][j].ne_j[1] = j-1;
      cells[i][j].ne_j[2] = j;
      cells[i][j].ne_j[3] = j+1;
      cells[i][j].ne_j[4] = j+1;
      cells[i][j].ne_j[5] = j+1;
      cells[i][j].ne_j[6] = j;
      cells[i][j].ne_j[7] = j-1;
      cells[i][j].ndist[0] = dx;
      cells[i][j].ndist[1] = dd;
      cells[i][j].ndist[2] = dy;
      cells[i][j].ndist[3] = dd;
      cells[i][j].ndist[4] = dx;
      cells[i][j].ndist[5] = dd;
      cells[i][j].ndist[6] = dy;
      cells[i][j].ndist[7] = dd;
    }
  }
}

void initiate_drainage_info_4(celltype **cells,meshtype mesh)
{

  int i,j;
  int ny = mesh.ny;
  int nx = mesh.nx;
  double dy = mesh.dy;
  double dx = mesh.dx;
  double dd = sqrt(dx*dx+dy*dy);

  /*corners*/

  /*lower left*/
  i = 0; j = 0;
  cells[i][j].nne = 2;
  cells[i][j].ne_i[0] = i+1;
  cells[i][j].ne_i[1] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j+1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dx;

  /*upper left*/
  i = ny-1; j = 0;
  cells[i][j].nne = 2;
  cells[i][j].ne_i[0] = i-1;
  cells[i][j].ne_i[1] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j+1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dx;

  /*upper right*/
  i = ny-1; j = nx-1;
  cells[i][j].nne = 2;
  cells[i][j].ne_i[0] = i-1;
  cells[i][j].ne_i[1] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j-1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dx;

  /*lower right*/
  i = 0; j = nx-1;
  cells[i][j].nne = 2;
  cells[i][j].ne_i[0] = i+1;
  cells[i][j].ne_i[1] = i;
  cells[i][j].ne_j[0] = j;
  cells[i][j].ne_j[1] = j-1;
  cells[i][j].ndist[0] = dy;
  cells[i][j].ndist[1] = dx;

  /*left side*/
  j = 0;
  for (i=1;i<ny-1;i++) {
    cells[i][j].nne = 3;
    cells[i][j].ne_i[0] = i-1;
    cells[i][j].ne_i[1] = i;
    cells[i][j].ne_i[2] = i+1;
    cells[i][j].ne_j[0] = j;
    cells[i][j].ne_j[1] = j+1;
    cells[i][j].ne_j[2] = j;
    cells[i][j].ndist[0] = dy;
    cells[i][j].ndist[1] = dx;
    cells[i][j].ndist[2] = dy;
  }

  /*right side*/
  j = nx-1;
  for (i=1;i<ny-1;i++) {
    cells[i][j].nne = 3;
    cells[i][j].ne_i[0] = i-1;
    cells[i][j].ne_i[1] = i;
    cells[i][j].ne_i[2] = i+1;
    cells[i][j].ne_j[0] = j;
    cells[i][j].ne_j[1] = j-1;
    cells[i][j].ne_j[2] = j;
    cells[i][j].ndist[0] = dy;
    cells[i][j].ndist[1] = dx;
    cells[i][j].ndist[2] = dy;
  }

  /*lower side*/
  i = 0;
  for (j=1;j<nx-1;j++) {
    cells[i][j].nne = 3;
    cells[i][j].ne_i[0] = i;
    cells[i][j].ne_i[1] = i+1;
    cells[i][j].ne_i[2] = i;
    cells[i][j].ne_j[0] = j-1;
    cells[i][j].ne_j[1] = j;
    cells[i][j].ne_j[2] = j+1;
    cells[i][j].ndist[0] = dx;
    cells[i][j].ndist[1] = dy;
    cells[i][j].ndist[2] = dx;
  }

  /*upper side*/
  i = ny-1;
  for (j=1;j<nx-1;j++) {
    cells[i][j].nne = 3;
    cells[i][j].ne_i[0] = i;
    cells[i][j].ne_i[1] = i-1;
    cells[i][j].ne_i[2] = i;
    cells[i][j].ne_j[0] = j-1;
    cells[i][j].ne_j[1] = j;
    cells[i][j].ne_j[2] = j+1;
    cells[i][j].ndist[0] = dx;
    cells[i][j].ndist[1] = dy;
    cells[i][j].ndist[2] = dx;
  }

  /*middle*/
  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {
      cells[i][j].nne = 4;
      cells[i][j].ne_i[0] = i+1;
      cells[i][j].ne_i[1] = i;
      cells[i][j].ne_i[2] = i-1;
      cells[i][j].ne_i[3] = i;
      cells[i][j].ne_j[0] = j;
      cells[i][j].ne_j[1] = j+1;
      cells[i][j].ne_j[2] = j;
      cells[i][j].ne_j[3] = j-1;
      cells[i][j].ndist[0] = dy;
      cells[i][j].ndist[1] = dx;
      cells[i][j].ndist[2] = dy;
      cells[i][j].ndist[3] = dx;
    }
  }
}


void fixpit(int ipit,int jpit,celltype **cells) 
{

  int k;
  int tn,ni,nj,ri,rj;
  double slope,mslope;
  double minz = 0.1;


  if (cells[ipit][jpit].fixflag <= 0) {

    /*find spillpoint neighbour*/
    tn = 0; ni = cells[ipit][jpit].ne_i[0]; nj = cells[ipit][jpit].ne_j[0];
    ri = ni; rj = nj;
    mslope = (cells[ipit][jpit].wsurf-cells[ni][nj].wsurf)/cells[ipit][jpit].ndist[0];
    for (k=1;k<cells[ipit][jpit].nne;k++) {
      ni = cells[ipit][jpit].ne_i[k]; nj = cells[ipit][jpit].ne_j[k];
      slope = (cells[ipit][jpit].wsurf-cells[ni][nj].wsurf)/cells[ipit][jpit].ndist[k];
      if (slope > mslope) {
	tn = k;
	ri = ni;
	rj = nj;
	mslope = slope;
      }/*if*/
    }/*k*/
  
    /*if at boundary*/
    /*    if (cells[ri][rj].fixflag > 0) {
      
      cells[ipit][jpit].reciever = -1;
      
      }*//*if*/
    
    /*if pit*/
    if (mslope <= 0.0) {
      
      /*elevate cell*/
      cells[ipit][jpit].wsurf = cells[cells[ipit][jpit].ne_i[tn]][cells[ipit][jpit].ne_j[tn]].wsurf + minz;
      
      /*registre reciever*/
      cells[ipit][jpit].reciever = tn;
      cells[ipit][jpit].ri = cells[ipit][jpit].ne_i[tn];
      cells[ipit][jpit].rj = cells[ipit][jpit].ne_j[tn];
      cells[ipit][jpit].rdist = cells[ipit][jpit].ndist[tn];

      
      /*call fixpit recursively*/
      for (k=0;k<cells[ipit][jpit].nne;k++) { 
	fixpit(cells[ipit][jpit].ne_i[k],cells[ipit][jpit].ne_j[k],cells);
      }
      
    }/*if*/
    
    else {
      
      /*registre reciever*/
      cells[ipit][jpit].reciever = tn;
      cells[ipit][jpit].ri = cells[ipit][jpit].ne_i[tn];
      cells[ipit][jpit].rj = cells[ipit][jpit].ne_j[tn];
      cells[ipit][jpit].rdist = cells[ipit][jpit].ndist[tn];
      
    }/*else*/
    
  }

}


void get_dnetwork(celltype **cells,meshtype *mesh)
{

  long i,j,k;
  long tn,ni,nj,ri,rj;
  double slope,mslope,mdist;
  int npit,*pit_i,*pit_j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;

  pit_i = (int*) malloc(nx*ny*sizeof(int));
  pit_j = (int*) malloc(nx*ny*sizeof(int));

  /*initialize*/
  npit = 0;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].reciever = -1;
      cells[i][j].wsurf = cells[i][j].bed;
      /*cells[i][j].wsurf = cells[i][j].topsedi+cells[i][j].Pw;*/
    }/*j*/
  }/*i*/

  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      /*if in drainage area*/
      if (cells[i][j].fixflag <= 0) {
    
	/*find steepest neighbour*/
	tn = 0; ni = cells[i][j].ne_i[0]; nj = cells[i][j].ne_j[0];
	ri = ni; rj = nj;
	mslope = (cells[i][j].wsurf-cells[ni][nj].wsurf)/cells[i][j].ndist[0];
	for (k=1;k<cells[i][j].nne;k++) {
	  ni = cells[i][j].ne_i[k]; nj = cells[i][j].ne_j[k];
	  slope = (cells[i][j].wsurf-cells[ni][nj].wsurf)/cells[i][j].ndist[k];
	  if (slope > mslope) {
	    tn = k;
	    ri = ni;
	    rj = nj;
	    mdist = cells[i][j].ndist[k];
	    mslope = slope;
	  }/*if*/
	}/*k*/

	if (mslope <= 0.0) {
	  
	  /*registre pit*/
	  pit_i[npit] = i;
	  pit_j[npit] = j;
	  npit += 1;

	}/*else*/
	
	/*if not pit*/
	else {
	  
	  /*registre reciever*/
	  cells[i][j].reciever = tn;
	  cells[i][j].ri = cells[i][j].ne_i[tn];
	  cells[i][j].rj = cells[i][j].ne_j[tn];
	  cells[i][j].rdist = cells[i][j].ndist[tn];
	  
	}/*else*/

      }/*if*/

    }/*j*/      
    
  }/*i*/
  
  /*printf("npit = %d\n",npit);*/


  /*if pits exist*/
  for (i=0;i<npit;i++) {
    
    fixpit(pit_i[i],pit_j[i],cells);
    
  }/*i*/


  /*compute lakewater*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].lakewater = cells[i][j].wsurf - cells[i][j].topice;
    }/*j*/
  }/*i*/


  free(pit_i);
  free(pit_j);

}/*get_dnetwork*/


void get_cascade(celltype **cells,meshtype mesh,long **cascade,long **catchment,long *ncascade,long *ncatchment) 
{

  long i,j,di,dj,ri,rj;
  long nstack,newstack,newcascade;
  int nx = mesh.nx;
  int ny = mesh.ny;
  long nc = (long)nx*(long)ny;
  int **parcels;
  int **stack;

  parcels =  (int **) malloc(ny*sizeof(int*)); for (i=0;i<ny;i++) parcels[i] = (int*) malloc(nx*sizeof(int));
  stack = (int **) malloc(nc*sizeof(int*)); for (i=0;i<nc;i++) stack[i] = (int*) malloc(2*sizeof(int));

  /*initiate*/
  (*ncascade) = 0;
  nstack = 0;
  newcascade = 0;
  newstack = 0;

  /*********** the first go ***************/

  /*give parcels*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if (cells[i][j].fixflag <= 0) parcels[i][j] = 1;
      else parcels[i][j] = 0;
    }
  }

  /*pass parcels*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if ((cells[i][j].reciever > -1)&&(cells[i][j].fixflag <= 0)) {
	parcels[cells[i][j].ri][cells[i][j].rj] += 1;
	parcels[i][j] -= 1;
      }
    }
  }

  /*place in cascade or stack*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      if (cells[i][j].fixflag <= 0) {

	if (parcels[i][j] == 0) {
	  (*ncascade) += 1;
	  cascade[(*ncascade)-1][0] = i;
	  cascade[(*ncascade)-1][1] = j;
	}
	else {
	  nstack += 1;
	  stack[nstack-1][0] = i;
	  stack[nstack-1][1] = j;
	}
      }
      
    }
  }

  while (nstack > 0) {

    /*initalize*/
    newstack = 0;
    newcascade = 0;

    /*give parcels*/
    for (i=0;i<nstack;i++) parcels[stack[i][0]][stack[i][1]] = 1;

    /*pass parcels*/
    for (i=0;i<nstack;i++) {
      di = stack[i][0]; dj = stack[i][1];
      if ((cells[di][dj].reciever > -1)&&(cells[di][dj].fixflag <= 0)) {
	ri = cells[di][dj].ri;
	rj = cells[di][dj].rj;
	parcels[ri][rj] += 1;
	parcels[di][dj] -= 1;
      }
    }

    /*place in cascade or stack*/
    for (i=0;i<nstack;i++) {

      di = stack[i][0];
      dj = stack[i][1];

      if (parcels[di][dj] == 0) {
	(*ncascade) += 1;
	newcascade += 1;
	cascade[(*ncascade)-1][0] = di;
	cascade[(*ncascade)-1][1] = dj;
      }
      else {
	newstack += 1;
	stack[newstack-1][0] = di;
	stack[newstack-1][1] = dj;
      }

    }

    /*if only pits left*/
    if (newcascade == 0) {
      for (i=0;i<newstack;i++) {
	(*ncascade) += 1;
	cascade[(*ncascade)-1][0] = stack[i][0];
	cascade[(*ncascade)-1][1] = stack[i][1];
      }
      nstack = 0;
    }
    else nstack = newstack;
    
  }

  /******* compute stream orders and parents ***********/

  /*initialize*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].nparent = 0;
      cells[i][j].streamorder = 1;
    }
  }
  (*ncatchment) = 0;

  /*loop through cascade*/
  for (i=0;i<(*ncascade);i++) {

    /*doner index*/ 
    di = cascade[i][0];
    dj = cascade[i][1];

    if (cells[di][dj].reciever > -1) {

      /*reciever index*/
      ri = cells[di][dj].ri;
      rj = cells[di][dj].rj;

      /*add stream order*/
      cells[ri][rj].streamorder += cells[di][dj].streamorder;
      cells[ri][rj].nparent += 1;
      cells[ri][rj].pa_i[cells[ri][rj].nparent-1] = di;
      cells[ri][rj].pa_j[cells[ri][rj].nparent-1] = dj;

    }/*if*/

    else {

      (*ncatchment) += 1;
      catchment[(*ncatchment)-1][0] = di;
      catchment[(*ncatchment)-1][1] = dj;

    }/*else*/

  }/*i*/


  for (i=0;i<ny;i++) free(parcels[i]); free(parcels);
  for (i=0;i<nc;i++) free(stack[i]); free(stack);

}


void add_to_catchment(int di,int dj,celltype **cells,int cnr)
{

  int i;

  /*registre catchment number of this cell*/
  cells[di][dj].cnumber = cnr;

  /*add parents to catchment recursively*/
  for (i=0;i<cells[di][dj].nparent;i++) add_to_catchment(cells[di][dj].pa_i[i],cells[di][dj].pa_j[i],cells,cnr);

}


void get_catchment(long ncatchment,long **catchment,celltype **cells)
{

  int di,dj;
  int cnr;

  /*loop catchment endpoints*/
  for (cnr=1;cnr<=ncatchment;cnr++) {

    /*start cell*/
    di = catchment[cnr-1][0];
    dj = catchment[cnr-1][1];

    /*add to catchment*/
    add_to_catchment(di,dj,cells,cnr);

  }/*cnr*/

}

void pass_water(celltype **cells,meshtype *mesh,fproptype fprop,long ncascade,long **cascade,double dt)
{

  int i,j;
  int di,dj,ri,rj,pi,pj;
  long k,tp;
  double Qw;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double area = (*mesh).dx*(*mesh).dy;
  double spy = 3600.0*24.0*365.25;

  /*initialize*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].water = cells[i][j].water + 1.0*dt;
      cells[i][j].Qw = 0.0;
    }
  }

  /*loop cascade and pass water*/
  for (k=0;k<ncascade;k++) {

    /*doner index*/
    di = cascade[k][0];
    dj = cascade[k][1];

    if (cells[di][dj].reciever >= 0) {

      /*reciever index*/
      ri = cells[di][dj].ri;
      rj = cells[di][dj].rj;
      
      /*water flux  m3/s */
      Qw = area*cells[di][dj].water/(spy*dt);
            
      /*store variables*/
      cells[di][dj].Qw = Qw;
      
      /*transfer water*/
      cells[ri][rj].water += cells[di][dj].water;


    }/*if*/

  }/*k*/

}

void smooth_water_flux(celltype **cells,meshtype mesh,cornertype **cp)
{

  int i,j;
  int nx = mesh.nx;
  int ny = mesh.ny;

#pragma omp parallel shared(cells,cp) private(i,j) firstprivate(nx,ny)
  {

#pragma omp for schedule(static)
  for (i=1;i<ny;i++) {
    for (j=1;j<nx;j++) {
      cp[i][j].Qw = 0.25*(cells[i-1][j-1].Qw+cells[i][j-1].Qw+cells[i][j].Qw+cells[i-1][j].Qw);
    }
  }

#pragma omp for schedule(static)
  for (i=1;i<(ny-1);i++) {
    for (j=1;j<(nx-1);j++) {
      cells[i][j].Qw = 0.25*(cp[i][j].Qw+cp[i+1][j].Qw+cp[i+1][j+1].Qw+cp[i][j+1].Qw);
    }
  }

  }


}


void fluvial_transport(celltype **cells,meshtype *mesh,fproptype fprop,long ncascade,long **cascade,double dt)
{

  int i,j;
  int di,dj,ri,rj,pi,pj;
  long k,tp;
  double Qw,slope,Wc,tau,Qt,sediload;
  double depo,ero,minptopo,maxdepo,coverfac;
  double mean_dHs = 0.0;
  double mean_erate = 0.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double area = (*mesh).dx*(*mesh).dy;

  double pr = fprop.pr;
  double rho_s = fprop.rho_s;
  double Dg = fprop.Dg;
  double kw = fprop.kw;
  double tau_c = fprop.tau_c;
  double Kt = fprop.Kt;
  double Ke = fprop.Ke;

  double rho_w = 1000.0;
  double spy = 3600.0*24.0*365.25;
  double g = 10.0;
  double Rb = (rho_s-rho_w)/rho_w;

  double cw = 1.0/2.0;
  double sw = 1.0/4.0;

  /*initialize*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].water = pr*dt;
      cells[i][j].Qw = 0.0;
      cells[i][j].Wc = 0.0;
      cells[i][j].tau = 0.0;
      cells[i][j].Qt = 0.0;
      cells[i][j].streamload = 0.0;
      cells[i][j].fluvial_erate = 0.0;
    }
  }

  /*loop cascade and compute carrying capacity*/
  for (k=0;k<ncascade;k++) {


    /*doner index*/
    di = cascade[k][0];
    dj = cascade[k][1];

    if (cells[di][dj].reciever >= 0) {

      /*reciever index*/
      ri = cells[di][dj].ri;
      rj = cells[di][dj].rj;
      
      /*water flux  m3/s */
      Qw = area*cells[di][dj].water/(spy*dt);
      
      /*channel width m*/
      Wc = kw*sqrt(Qw);
      
      /*bed slope*/
      slope = (cells[di][dj].topsedi-cells[ri][rj].topsedi)/cells[di][dj].rdist;      

      /*if downhill*/
      if (slope > 0.0) {
	
	tau = rho_w*pow(g*Qw*slope,2.0/3.0)*pow(Wc,-2.0/3.0);
	tau /= (rho_s-rho_w)*g*Dg;
	
	if (tau > tau_c) Qt = Kt*Wc/(Rb*g)*pow(tau/rho_w,1.5)*pow(1.0-0.846*tau_c/tau,4.5);
	else Qt = 0.0;
	
      }
      else {
	tau = 0.0;
	Qt = 0.0;
      }
      
      /*store variables*/
      cells[di][dj].Qw = Qw;
      cells[di][dj].Wc = Wc;
      cells[di][dj].tau = tau;
      cells[di][dj].Qt = Qt;
      
      /*transfer water*/
      cells[ri][rj].water += cells[di][dj].water;

      /*transfer sediment*/
      cells[ri][rj].streamload += cells[di][dj].streamload;
      mean_dHs += cells[di][dj].streamload;
      cells[di][dj].streamload = 0.0;
      

    }/*if*/

  }/*k*/

  /*loop cascade again and transport sediment*/
  for (k=0;k<ncascade;k++) {


    /*doner index*/
    di = cascade[k][0];
    dj = cascade[k][1];

    if (cells[di][dj].reciever >= 0) {

      /*reciever index*/
      ri = cells[di][dj].ri;
      rj = cells[di][dj].rj;
      
      Qt = cells[di][dj].Qt;
      if (cells[ri][rj].Qt < Qt) Qt = cells[ri][rj].Qt;
	
      /*sediment thickness capacity*/
      sediload = Qt*spy*dt/area;

	
      /*entrain*/
      if ((sediload > cells[di][dj].streamload)&&(sediload > 0.0)) {

	ero = sediload - cells[di][dj].streamload;
	if (ero > cells[di][dj].sedi) ero = cells[di][dj].sedi;
	if (ero > (cells[di][dj].topsedi-cells[ri][rj].bed)) ero = cells[di][dj].topsedi-cells[ri][rj].bed;

	if (ero < 0.0) ero = 0.0;
	cells[di][dj].sedi -= ero;
	cells[di][dj].streamload += ero;


	/*consider bedrock erosion*/

	/*water flux  m3/s */
	Qw = area*cells[di][dj].water/(spy*dt);
      
	/*channel width m*/
	Wc = kw*sqrt(Qw);
      
	/*bed slope*/
	slope = (cells[di][dj].bed-cells[ri][rj].bed)/cells[di][dj].rdist;      

	/*if downhill*/
	if (slope > 0.0) {
	
	  tau = rho_w*pow(g*Qw*slope,2.0/3.0)*pow(Wc,-2.0/3.0);
	  tau /= (rho_s-rho_w)*g*Dg;
	
	  /*bedrock erosion*/
	  if ((sediload > cells[di][dj].streamload)&&(tau > tau_c)) {

	    coverfac = (sediload - cells[di][dj].streamload)/sediload;
	    /*ero = Ke*coverfac*pow(Qw,.5)*slope;*/
	    ero = Ke*coverfac*(tau-tau_c);
	    
	    if ((di > 0)&&(dj > 0)&&(di < ny-1)&&(dj < nx-1)) {
	      
	      if ((ri == di+1)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di+1)&&(rj == dj)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di+1)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		mean_erate += cw*ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj].streamload += cw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
	      }
	    }
	    else {
	      if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
	      if (ero < 0.0) ero = 0.0;
	      cells[di][dj].fluvial_erate = cw*ero;
	      mean_erate += cw*ero;
	      ero *= dt;
	      cells[di][dj].bed -= cw*ero;
	      cells[di][dj].streamload += cw*ero;
	      cells[di][dj].fluvial_erosion += cw*ero;
	    }
	    
	  }/*if erosion*/
	  
	}/*if down slope*/
      }/*if*/
	
      /*deposit*/
      else {
	
	depo = cells[di][dj].streamload - sediload;

	/*min parent topography*/
	if (cells[di][dj].nparent > 0) {
	  pi = cells[di][dj].pa_i[0];
	  pj = cells[di][dj].pa_j[0];
	  minptopo = cells[pi][pj].topsedi;
	  for (i=1;i<cells[di][dj].nparent;i++) {
	    pi = cells[di][dj].pa_i[i];
	    pj = cells[di][dj].pa_j[i];
	    if (minptopo > cells[pi][pj].topsedi) minptopo = cells[pi][pj].topsedi;
	  }
	}
	maxdepo = minptopo - cells[di][dj].topsedi;
	if (maxdepo > 0.0) {
	  if (depo > maxdepo) depo = maxdepo;
	  cells[di][dj].sedi += depo;
	  cells[di][dj].streamload -= depo;
	}
      }

      /*transfer sediment*/
      cells[ri][rj].streamload += cells[di][dj].streamload;
      mean_dHs += cells[di][dj].streamload;
      cells[di][dj].streamload = 0.0;
      

    }/*if*/

  }/*k*/

  (*mesh).mean_dHs_fluvial = mean_dHs/((double)(*mesh).nc*dt);
  (*mesh).mean_fluvial_erate = mean_erate/((double)(*mesh).nc);

}

void fluvial_transport2(celltype **cells,meshtype *mesh,fproptype fprop,long ncascade,long **cascade,double dt)
{

  int i,j;
  int di,dj,ri,rj,pi,pj;
  long k,tp;
  double Qw,slope,Wc,tau,Qt,sediload;
  double depo,ero,minptopo,maxdepo,coverfac;
  double mean_dHs = 0.0;
  double mean_erate = 0.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double area = (*mesh).dx*(*mesh).dy;

  double pr = fprop.pr;
  double rho_s = fprop.rho_s;
  double Dg = fprop.Dg;
  double kw = fprop.kw;
  double tau_c = fprop.tau_c;
  double Kt = fprop.Kt;
  double Ke = fprop.Ke;

  double rho_w = 1000.0;
  double spy = 3600.0*24.0*365.25;
  double g = 10.0;
  double Rb = (rho_s-rho_w)/rho_w;

  double cw = 1.0/2.0;
  double sw = 1.0/4.0;

  /*initialize*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].water = pr*dt;
      cells[i][j].Qw = 0.0;
      cells[i][j].Wc = 0.0;
      cells[i][j].tau = 0.0;
      cells[i][j].Qt = 0.0;
      cells[i][j].streamload = 0.0;
      cells[i][j].fluvial_erate = 0.0;
    }
  }

  /*loop cascade and compute carrying capacity*/
  for (k=0;k<ncascade;k++) {


    /*doner index*/
    di = cascade[k][0];
    dj = cascade[k][1];

    if (cells[di][dj].reciever >= 0) {

      /*reciever index*/
      ri = cells[di][dj].ri;
      rj = cells[di][dj].rj;
      
      /*water flux  m3/s */
      Qw = area*cells[di][dj].water/(spy*dt);
      
      /*channel width m*/
      Wc = kw*sqrt(Qw);
      
      /*bed slope*/
      slope = (cells[di][dj].topsedi-cells[ri][rj].topsedi)/cells[di][dj].rdist;      


      if (slope > 0.0) {
	
	tau = rho_w*pow(g*Qw*slope,2.0/3.0)*pow(Wc,-2.0/3.0);
	tau /= (rho_s-rho_w)*g*Dg;
	
	if (tau > tau_c) Qt = Kt*Wc/(Rb*g)*pow(tau/rho_w,1.5)*pow(1.0-0.846*tau_c/tau,4.5);
	else Qt = 0.0;
	
	/*sediment thickness capacity*/
	sediload = Qt*spy*dt/area;
	
	/*entrain*/
	if ((sediload > cells[di][dj].streamload)&&(sediload > 0.0)) {

	  ero = sediload - cells[di][dj].streamload;
	  if (ero > cells[di][dj].sedi) ero = cells[di][dj].sedi;
	  if (ero < 0.0) ero = 0.0;
	  cells[di][dj].sedi -= ero;
	  cells[di][dj].streamload += ero;

	  /*bedrock erosion*/
	  if (sediload > cells[di][dj].streamload) {

	    coverfac = (sediload - cells[di][dj].streamload)/sediload;
	    ero = Ke*coverfac*pow(Qw,.5)*slope;
	  
	    if ((di > 0)&&(dj > 0)&&(di < ny-1)&&(dj < nx-1)) {

	      if ((ri == di+1)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di+1)&&(rj == dj)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di+1)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj+1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj+1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		cells[di][dj+1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj+1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj+1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di-1)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di][dj-1].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di][dj-1].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di][dj-1].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di][dj-1].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else if ((ri == di)&&(rj == dj-1)) {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero > (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di-1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero > (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt)) ero = (cells[di+1][dj].bed-cells[ri][rj].bed)/(sw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		cells[di-1][dj].fluvial_erate = sw*ero;
		cells[di+1][dj].fluvial_erate = sw*ero;
		mean_erate += ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di-1][dj].bed -= sw*ero;
		cells[di+1][dj].bed -= sw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
		cells[di-1][dj].fluvial_erosion += sw*ero;
		cells[di+1][dj].fluvial_erosion += sw*ero;
		cells[di][dj].streamload += ero;
	      }
	      else {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		mean_erate += cw*ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj].streamload += cw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
	      }
	    }
	    else {
		if (ero > (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt)) ero = (cells[di][dj].bed-cells[ri][rj].bed)/(cw*dt);
		if (ero < 0.0) ero = 0.0;
		cells[di][dj].fluvial_erate = cw*ero;
		mean_erate += cw*ero;
		ero *= dt;
		cells[di][dj].bed -= cw*ero;
		cells[di][dj].streamload += cw*ero;
		cells[di][dj].fluvial_erosion += cw*ero;
	    }

	  }/*if*/
	}/*if*/
	/*deposit*/
	else {

	  depo = cells[di][dj].streamload - sediload;

	  /*min parent topography*/
	  if (cells[di][dj].nparent > 0) {
	    pi = cells[di][dj].pa_i[0];
	    pj = cells[di][dj].pa_j[0];
	    minptopo = cells[pi][pj].topsedi;
	    for (i=1;i<cells[di][dj].nparent;i++) {
	      pi = cells[di][dj].pa_i[i];
	      pj = cells[di][dj].pa_j[i];
	      if (minptopo > cells[pi][pj].topsedi) minptopo = cells[pi][pj].topsedi;
	    }
	  }
	  maxdepo = minptopo - cells[di][dj].topsedi;
	  if (maxdepo > 0.0) {
	    if (depo > maxdepo) depo = maxdepo;
	    cells[di][dj].sedi += depo;
	    cells[di][dj].streamload -= depo;
	    }
	}
      }
      else {
	tau = 0.0;
	Qt = 0.0;
      }
      
      /*store variables*/
      cells[di][dj].Qw = Qw;
      cells[di][dj].Wc = Wc;
      cells[di][dj].tau = tau;
      cells[di][dj].Qt = Qt;
      
      /*transfer water*/
      cells[ri][rj].water += cells[di][dj].water;

      /*transfer sediment*/
      cells[ri][rj].streamload += cells[di][dj].streamload;
      mean_dHs += cells[di][dj].streamload;
      cells[di][dj].streamload = 0.0;
      

    }/*if*/

  }/*k*/

  (*mesh).mean_dHs_fluvial = mean_dHs/((double)(*mesh).nc*dt);
  (*mesh).mean_fluvial_erate = mean_erate/((double)(*mesh).nc);

}
