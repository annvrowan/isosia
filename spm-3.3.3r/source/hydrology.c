

double max(double val1,double val2)
{

  double maxval = val1;
  if (val2 > maxval) maxval = val2;

  return(maxval);

}


void glacial_hydrology(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hwproptype hwprop,double dt)
{

  int i,j;
  double dHw,Hw,dphi;
  double meanpw = 0.0;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double po = hwprop.po; /*ice porosity*/

#pragma omp parallel shared(cells,hp,vp,meanpw) private(i,j,Hw,dHw,dphi) firstprivate(nx,ny,dx,dy,po)
  {

    /*compute pressure and initiate dHw*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	/*penalties*/
	if ((cells[i][j].ice < 1.0)&&(cells[i][j].Hgw > 0.0)) cells[i][j].Hgw = 0.0;
	if ((cells[i][j].Hgw > po*cells[i][j].tn)&&(cells[i][j].tn > 0.0)) cells[i][j].Hgw = po*cells[i][j].tn;
	/*water pressure function*/
	cells[i][j].Pw = cells[i][j].Hgw/po;
	cells[i][j].dHw = 0.0; /*initialize*/
      }
    }

    /*loop h-points*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	Hw = max(cells[i][j-1].Hgw,cells[i][j].Hgw); if (Hw < 0.0) Hw = 0.0;
	dphi = hp[i][j].dbdx + (cells[i][j].Pw-cells[i][j-1].Pw)/dx;
	dHw = -hp[i][j].Kgw*Hw*dphi*dt/dx;
	if (dHw < 0.0) {
	  if (cells[i][j].Hgw <= 0.0) dHw = 0.0;
	  else if (dHw < -cells[i][j].Hgw) dHw = -cells[i][j].Hgw;
	}
	else {
	  if (cells[i][j-1].Hgw <= 0.0) dHw = 0.0;
	  else if (dHw > cells[i][j-1].Hgw) dHw = cells[i][j-1].Hgw;
	}
	hp[i][j].dHw = dHw;
	/*#pragma omp critical*/
	/*{
	  cells[i][j].dHw += dHw;
	  cells[i][j-1].dHw -= dHw;
	  }*/
      }
    }

    /*loop v-points*/
#pragma omp for schedule(static)
    for (j=0;j<nx;j++) {
      for (i=1;i<ny;i++) {
	Hw = max(cells[i-1][j].Hgw,cells[i][j].Hgw); if (Hw <= 0.0) Hw = 0.0;
	dphi = vp[i][j].dbdy + (cells[i][j].Pw-cells[i-1][j].Pw)/dy; 
	dHw = -vp[i][j].Kgw*Hw*dphi*dt/dy;
	if (dHw < 0.0) {
	  if (cells[i][j].Hgw <= 0.0) dHw = 0.0; 
	  else if (dHw < -cells[i][j].Hgw) dHw = -cells[i][j].Hgw;
	}
	else {
	  if (cells[i-1][j].Hgw <= 0.0) dHw = 0.0;
	  else if (dHw > cells[i-1][j].Hgw) dHw = cells[i-1][j].Hgw;
	}
	vp[i][j].dHw = dHw;
	/*	#pragma omp critical*/
	/*{
	  cells[i][j].dHw += dHw;
	  cells[i-1][j].dHw -= dHw;
	  }*/
      }
    }

#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	/*cells[i][j].Hgw = cells[i][j].Hgw + cells[i][j].dHw;*/
	cells[i][j].Hgw = cells[i][j].Hgw + hp[i][j].dHw - hp[i][j+1].dHw + vp[i][j].dHw - vp[i+1][j].dHw;
      }
    }

    /*compute mean water pressure*/
#pragma omp for schedule(static) reduction(+:meanpw) 
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	meanpw += cells[i][j].Hgw/po;
      }
    }

#pragma omp single
    {
      meanpw /= (double)((*mesh).nc);
      (*mesh).meanpw = meanpw;
    }    


  }/*pragma*/

}/*void*/



void glacial_hydrology3(celltype **cells,meshtype *mesh,hwproptype hwprop)
{

  int i,j,step;
  double pval,dHw,Hw,dphi,fac1,fac2;
  double res,totHgw,residual,meanpw;
  double h1,h2,h3,h4;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double Kgw = hwprop.Kgw;
  double po = hwprop.po; /*ice porosity*/

  double minPw = 1.0e-2;
  double mindH = 1.0e-0;

  
  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {
      cells[i][j].Hgw = cells[i][j].bed + cells[i][j].Pw;
    }
  }



  fac1 = 2.0/(dx*dx)+2.0/(dy*dy); 
  residual = 1.0; step = 1;

  while ((residual > 1.0e-3)&&(step < 1000)) {

    res = 0.0; totHgw = 0.0;
    for (i=1;i<ny-1;i++) {
      for (j=1;j<nx-1;j++) {
	
	if (cells[i][j].margin == -1) {

	  fac2 = Kgw*(cells[i][j].Hgw-cells[i][j].bed + 1.0);

	  /*
	  if ((cells[i+1][j].Pw > minPw)||(cells[i+1][j].Hgw < cells[i][j].Hgw)) h1 = cells[i+1][j].Hgw; else h1 = cells[i][j].Hgw;
	  if ((cells[i][j+1].Pw > minPw)||(cells[i][j+1].Hgw < cells[i][j].Hgw)) h2 = cells[i][j+1].Hgw; else h2 = cells[i][j].Hgw;
	  if ((cells[i-1][j].Pw > minPw)||(cells[i-1][j].Hgw < cells[i][j].Hgw)) h3 = cells[i-1][j].Hgw; else h3 = cells[i][j].Hgw;
	  if ((cells[i][j-1].Pw > minPw)||(cells[i][j-1].Hgw < cells[i][j].Hgw)) h4 = cells[i][j-1].Hgw; else h4 = cells[i][j].Hgw;
	  */
	  if ((cells[i+1][j].dHw > mindH)||(cells[i+1][j].Hgw < cells[i][j].Hgw)) h1 = cells[i+1][j].Hgw; else h1 = cells[i][j].bed;
	  if ((cells[i][j+1].dHw > mindH)||(cells[i][j+1].Hgw < cells[i][j].Hgw)) h2 = cells[i][j+1].Hgw; else h2 = cells[i][j].bed;
	  if ((cells[i-1][j].dHw > mindH)||(cells[i-1][j].Hgw < cells[i][j].Hgw)) h3 = cells[i-1][j].Hgw; else h3 = cells[i][j].bed;
	  if ((cells[i][j-1].dHw > mindH)||(cells[i][j-1].Hgw < cells[i][j].Hgw)) h4 = cells[i][j-1].Hgw; else h4 = cells[i][j].bed;


	  pval = cells[i][j].dHw/(po*fac2)+(h4+h2)/(dx*dx)+(h3+h1)/(dy*dy)+((cells[i][j+1].Pw-cells[i][j-1].Pw)*(h2-h4)/(4.0*dx*dx)+(cells[i+1][j].Pw-cells[i-1][j].Pw)*(h1-h2)/(4.0*dy*dy))/fac2;
	  pval /= fac1;
	  if (pval < cells[i][j].bed) pval = cells[i][j].bed;
	  if (pval > cells[i][j].topice+100.0) pval = cells[i][j].topice+100.0;
	  res += pow(pval-cells[i][j].Hgw,2.0);
	  totHgw += pow(pval,2.0);
	  cells[i][j].Hgw = pval;
	  cells[i][j].Pw = pval - cells[i][j].bed;

	}
	else {
	  cells[i][j].Pw = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed-1.0;
	}

      }
    }

    for (i=0;i<ny;i++) 
      cells[i][0].Pw = cells[i][1].Pw;

    for (j=0;j<nx;j++) {
      cells[0][j].Pw = cells[1][j].Pw;
      cells[ny-1][j].Pw = cells[ny-2][j].Pw;
    }

    residual = sqrt(res/(totHgw+1.0e-6));

    /*printf("residual = %4.4e\n",residual); */
    step += 1;

  }

  meanpw = 0.0;
  /*compute mean water pressure*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanpw += cells[i][j].Pw;
    }
  }
  meanpw /= (double)((*mesh).nc);
  (*mesh).meanpw = meanpw;
  (*mesh).hydroitt = step;



}/*void*/

void glacial_hydrology4(celltype **cells,meshtype *mesh,hwproptype hwprop)
{

  int i,j,step;
  double pval,dHw,Hw,dphi,fac1,fac2;
  double res,totHgw,residual,meanpw;
  double h1,h2,h3,h4;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double Kgw = hwprop.Kgw;
  double po = hwprop.po; /*ice porosity*/

  double minPw = 1.0e-2;
  double minice = 20.0;
  double mindH = 1.0e-6;
  
  for (i=1;i<ny-1;i++) {
    for (j=1;j<nx-1;j++) {
      cells[i][j].Hgw = cells[i][j].bed + cells[i][j].Pw;
    }
  }

  fac1 = 2.0/(dx*dx)+2.0/(dy*dy); 
  residual = 1.0; step = 1;


  while ((residual > 1.0e-4)&&(step < 1000)) {

    res = 0.0; totHgw = 0.0;
    for (i=1;i<ny-1;i++) {
      for (j=1;j<nx-1;j++) {
	
	if (cells[i][j].margin == -1) {

	  fac2 = Kgw;

	  /*
	  if ((cells[i+1][j].Pw > minPw)||(cells[i+1][j].Hgw < cells[i][j].Hgw)) h1 = cells[i+1][j].Hgw; else h1 = cells[i][j].Hgw;
	  if ((cells[i][j+1].Pw > minPw)||(cells[i][j+1].Hgw < cells[i][j].Hgw)) h2 = cells[i][j+1].Hgw; else h2 = cells[i][j].Hgw;
	  if ((cells[i-1][j].Pw > minPw)||(cells[i-1][j].Hgw < cells[i][j].Hgw)) h3 = cells[i-1][j].Hgw; else h3 = cells[i][j].Hgw;
	  if ((cells[i][j-1].Pw > minPw)||(cells[i][j-1].Hgw < cells[i][j].Hgw)) h4 = cells[i][j-1].Hgw; else h4 = cells[i][j].Hgw;
	  */
	  
	  if ((cells[i+1][j].dHw > mindH)||(cells[i+1][j].Hgw < cells[i][j].Hgw)) h1 = cells[i+1][j].Hgw; else h1 = cells[i][j].bed;
	  if ((cells[i][j+1].dHw > mindH)||(cells[i][j+1].Hgw < cells[i][j].Hgw)) h2 = cells[i][j+1].Hgw; else h2 = cells[i][j].bed;
	  if ((cells[i-1][j].dHw > mindH)||(cells[i-1][j].Hgw < cells[i][j].Hgw)) h3 = cells[i-1][j].Hgw; else h3 = cells[i][j].bed;
	  if ((cells[i][j-1].dHw > mindH)||(cells[i][j-1].Hgw < cells[i][j].Hgw)) h4 = cells[i][j-1].Hgw; else h4 = cells[i][j].bed;
	  
	  /*
	  h1 = cells[i+1][j].Hgw;
	  h2 = cells[i][j+1].Hgw; 
	  h3 = cells[i-1][j].Hgw;
	  h4 = cells[i][j-1].Hgw;
	  */

	  pval = cells[i][j].dHw/fac2+(h4+h2)/(dx*dx)+(h3+h1)/(dy*dy);
	  pval /= fac1;
	  if (pval < cells[i][j].bed) pval = cells[i][j].bed;
	  if (pval > cells[i][j].topice+100.0) pval = cells[i][j].topice+100.0;
	  res += pow(pval-cells[i][j].Hgw,2.0);
	  totHgw += pow(pval,2.0);
	  cells[i][j].Hgw = pval;
	  cells[i][j].Pw = pval - cells[i][j].bed;

	}
	else {
	  cells[i][j].Pw = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed-1.0;
	}

      }
    }

    for (i=0;i<ny;i++) 
      cells[i][0].Pw = 0.5*cells[i][0].ice;
    /*cells[i][0].Pw = cells[i][1].Pw;*/

    for (j=0;j<nx;j++) {
      cells[0][j].Pw = cells[1][j].Pw;
      cells[ny-1][j].Pw = cells[ny-2][j].Pw;
    }

    residual = sqrt(res/(totHgw+1.0e-6));

    /*printf("residual = %4.4e\n",residual); */
    step += 1;

  }

  meanpw = 0.0;
  /*compute mean water pressure*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanpw += cells[i][j].Pw;
    }
  }
  meanpw /= (double)((*mesh).nc);
  (*mesh).meanpw = meanpw;
  (*mesh).hydroitt = step;

}/*void*/


