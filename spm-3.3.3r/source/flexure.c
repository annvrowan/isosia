
void isostasy(celltype **cells,meshtype *mesh,double **W)
{

  /*this solution converges slowly for large D values*/

  int i,j;
  int itt,Nitt = 20;
  double diff,diff_tot,W_tot,res,maxres;
  double mean_isostasy;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double scfac = 1.0e-6;
  double E = 70.0e9*scfac;
  double Te = 5000.0;
  double po = 0.25;
  double D = E*pow(Te,3.0)/(12.0*(1.0-po*po));
  double g = 10.0;
  double rho_a = 3200.0*scfac;
  double rho_r = 2900.0*scfac;
  double rho_s = 2500.0*scfac;
  double rho_i = 950.0*scfac;

  double fac = 6.0*D/pow(dx,4.0) + 8.0*D/(pow(dx,2.0)*pow(dy,2.0)) + 6.0*D/pow(dy,4.0) + rho_a*g;

  /*transfer existing isostatic deflections*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      W[i+2][j+2] = cells[i][j].isostasy;
    }
  }

  maxres = 1.0; itt = 0;
  while ((maxres > 1.0e-3)&&(itt < Nitt)) {

    maxres = 0.0;

    /*apply padding*/
    for (i=2;i<ny+2;i++) {
      W[i][0] = W[i][3];
      W[i][1] = W[i][2];
      W[i][nx+3] = W[i][nx];
      W[i][nx+2] = W[i][nx+1];
    }
    for (j=2;j<nx+2;j++) {
      W[0][j] = W[3][j];
      W[1][j] = W[2][j];
      W[ny+3][j] = W[ny][j];
      W[ny+2][j] = W[ny+1][j];
    }
    W[0][0] = W[3][3];
    W[0][1] = W[3][2];
    W[1][0] = W[2][3];
    W[1][1] = W[2][2];
    W[ny+2][0] = W[ny+1][3];
    W[ny+2][1] = W[ny+1][2];
    W[ny+3][0] = W[ny][3];
    W[ny+3][1] = W[ny][2];
    W[ny+2][nx+2] = W[ny+1][nx+1];
    W[ny+2][nx+3] = W[ny+1][nx];
    W[ny+3][nx+2] = W[ny][nx+1];
    W[ny+3][nx+3] = W[ny][nx];
    W[0][nx+2] = W[3][nx+1];
    W[0][nx+3] = W[3][nx];
    W[1][nx+2] = W[2][nx+1];
    W[1][nx+3] = W[2][nx];
    
    diff_tot = 0.0;
    W_tot = 0.0;
    for (i=2;i<ny+2;i++) {
      for (j=2;j<nx+2;j++) {
	
	diff = -W[i][j];

	W[i][j] = (rho_r*g*cells[i-2][j-2].abrasion-rho_i*g*cells[i-2][j-2].ice) + D*(4.0*(W[i][j+1]+W[i][j-1])-(W[i][j+2]+W[i][j-2]))/pow(dx,4.0) + 2.0*D*(2.0*(W[i][j+1]+W[i+1][j]+W[i][j-1]+W[i-1][j])-(W[i+1][j+1]+W[i+1][j-1]+W[i-1][j-1]+W[i-1][j+1]))/(pow(dx,2.0)*pow(dy,2.0)) + D*(4.0*(W[i+1][j]+W[i-1][j])-(W[i+2][j]+W[i-2][j]))/pow(dy,4.0);
	W[i][j] /= fac;

	diff += W[i][j];


	res = fabs(diff)/(fabs(W[i][j])+1.0e-16);
	if (res > maxres) maxres = res;

	diff_tot += diff*diff;
	W_tot += W[i][j]*W[i][j];

      }
    }


    /*apply padding*/
    for (i=2;i<ny+2;i++) {
      W[i][0] = W[i][3];
      W[i][1] = W[i][2];
      W[i][nx+3] = W[i][nx];
      W[i][nx+2] = W[i][nx+1];
    }
    for (j=2;j<nx+2;j++) {
      W[0][j] = W[3][j];
      W[1][j] = W[2][j];
      W[ny+3][j] = W[ny][j];
      W[ny+2][j] = W[ny+1][j];
    }
    W[0][0] = W[3][3];
    W[0][1] = W[3][2];
    W[1][0] = W[2][3];
    W[1][1] = W[2][2];
    W[ny+2][0] = W[ny+1][3];
    W[ny+2][1] = W[ny+1][2];
    W[ny+3][0] = W[ny][3];
    W[ny+3][1] = W[ny][2];
    W[ny+2][nx+2] = W[ny+1][nx+1];
    W[ny+2][nx+3] = W[ny+1][nx];
    W[ny+3][nx+2] = W[ny][nx+1];
    W[ny+3][nx+3] = W[ny][nx];
    W[0][nx+2] = W[3][nx+1];
    W[0][nx+3] = W[3][nx];
    W[1][nx+2] = W[2][nx+1];
    W[1][nx+3] = W[2][nx];
    
    for (i=ny+1;i>1;i--) {
      for (j=nx+1;j>1;j--) {
	
	diff = -W[i][j];

	W[i][j] = (rho_r*g*cells[i-2][j-2].abrasion-rho_i*g*cells[i-2][j-2].ice) + D*(4.0*(W[i][j+1]+W[i][j-1])-(W[i][j+2]+W[i][j-2]))/pow(dx,4.0) + 2.0*D*(2.0*(W[i][j+1]+W[i+1][j]+W[i][j-1]+W[i-1][j])-(W[i+1][j+1]+W[i+1][j-1]+W[i-1][j-1]+W[i-1][j+1]))/(pow(dx,2.0)*pow(dy,2.0)) + D*(4.0*(W[i+1][j]+W[i-1][j])-(W[i+2][j]+W[i-2][j]))/pow(dy,4.0);
	W[i][j] /= fac;

	diff += W[i][j];

	res = fabs(diff)/(fabs(W[i][j])+1.0e-16);
	if (res > maxres) maxres = res;

	diff_tot += diff*diff;
	W_tot += W[i][j]*W[i][j];

      }
    }


    maxres = sqrt(diff_tot/(W_tot + 1.0e-16));
    /*printf("itt = %d, maxres = %4.4e, W[50][50] = %4.4e\n",itt,maxres,W[50][50]);*/
    itt += 1;
    
  }


  /*save new isostasy*/
  mean_isostasy = 0.0;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      cells[i][j].isostasy = W[i+2][j+2];
      mean_isostasy += W[i+1][j+2];
    }
  }
  (*mesh).mean_isostasy = mean_isostasy/(ny*nx);

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
  double rho_r = 2600.0;
  double rho_s = 2500.0;
  double rho_i = 950.0;

  mean_isostasy = 0.0;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      
      load = -rho_i*g*cells[i][j].ice+rho_r*g*(cells[i][j].periglacial_erosion+cells[i][j].abrasion+cells[i][j].quarrying+cells[i][j].fluvial_erosion+cells[i][j].landslide_erosion+cells[i][j].hillslope_erosion);
      iso = load/(rho_a*g);				   
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

