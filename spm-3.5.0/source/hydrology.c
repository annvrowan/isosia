

double max(double val1,double val2)
{

  double maxval = val1;
  if (val2 > maxval) maxval = val2;

  return(maxval);

}


int compare_values(const void *elem1, const void *elem2)
{

  sortarraytype *i1,*i2;
  i1 = (sortarraytype*)elem1;
  i2 = (sortarraytype*)elem2;
  if (i1->value < i2->value) return 1; 
  else if (i1->value == i2->value) return 0;
  else return -1;

}


void get_discharge(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,double dt)
{

  int i,j; 
  double vout;



  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nn = nx*ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double spy = 365.25*24.0*3600.0;
  double minice = 10.0;
  double Tfac = 0.1;  

  /*Compute head*/
  for (i=0;i<ny;i++) { 
    for (j=0;j<nx;j++) {
      cells[i][j].Hgw = cells[i][j].topice;
    }
  }


  /*loop cells*/ 
  for (i=0;i<ny;i++) { 
    for (j=0;j<nx;j++) {
 
      if (cells[i][j].ice > minice) {
      
      cells[i][j].Qstack += cells[i][j].water*dt;

      /*vout = 0.0; 
      if (j > 0) if ((hp[i][j].vx < 0.0)&&(cells[i][j-1].ice > minice)) vout -= hp[i][j].vx;
      if (j < nx-1) if ((hp[i][j+1].vx > 0.0)&&(cells[i][j+1].ice > minice)) vout += hp[i][j+1].vx;
      if (i > 0) if ((vp[i][j].vy < 0.0)&&(cells[i-1][j].ice > minice)) vout -= vp[i][j].vy;
      if (i < ny-1) if ((vp[i+1][j].vy > 0.0)&&(cells[i+1][j].ice > minice)) vout += vp[i+1][j].vy;
     
      if (vout > 0.0) {
	if (j > 0) if ((hp[i][j].vx < 0.0)&&(cells[i][j-1].ice > minice)) cells[i][j-1].Qstack += cells[i][j].Qstack*(-hp[i][j].vx)/vout;
	if (j < nx-1) if ((hp[i][j+1].vx > 0.0)&&(cells[i][j+1].ice > minice)) cells[i][j+1].Qstack += cells[i][j].Qstack*(hp[i][j+1].vx)/vout;
	if (i > 0) if ((vp[i][j].vy < 0.0)&&(cells[i-1][j].ice > minice)) cells[i-1][j].Qstack += cells[i][j].Qstack*(-vp[i][j].vy)/vout;
	if (i < ny-1) if ((vp[i+1][j].vy > 0.0)&&(cells[i+1][j].ice > minice)) cells[i+1][j].Qstack += cells[i][j].Qstack*(vp[i+1][j].vy)/vout;
      } 
      else {*/

      vout = 0.0;
      if (j > 0) if (cells[i][j-1].Hgw < cells[i][j].Hgw) vout += cells[i][j-1].Kgw*(cells[i][j].Hgw - cells[i][j-1].Hgw);
      if (j < nx-1) if (cells[i][j+1].Hgw < cells[i][j].Hgw) vout += cells[i][j+1].Kgw*(cells[i][j].Hgw - cells[i][j+1].Hgw);
      if (i > 0) if (cells[i-1][j].Hgw < cells[i][j].Hgw) vout += cells[i-1][j].Kgw*(cells[i][j].Hgw - cells[i-1][j].Hgw);
      if (i < ny-1) if (cells[i+1][j].Hgw < cells[i][j].Hgw) vout += cells[i+1][j].Kgw*(cells[i][j].Hgw - cells[i+1][j].Hgw);
      
      if (vout > 0.0) {
	if (j > 0) if (cells[i][j-1].Hgw < cells[i][j].Hgw) cells[i][j-1].Qstack += cells[i][j].Qstack*cells[i][j-1].Kgw*(cells[i][j].Hgw-cells[i][j-1].Hgw)/vout;
	if (j < nx-1) if (cells[i][j+1].Hgw < cells[i][j].Hgw) cells[i][j+1].Qstack += cells[i][j].Qstack*cells[i][j+1].Kgw*(cells[i][j].Hgw-cells[i][j+1].Hgw)/vout;
	if (i > 0) if (cells[i-1][j].Hgw < cells[i][j].Hgw) cells[i-1][j].Qstack += cells[i][j].Qstack*cells[i-1][j].Kgw*(cells[i][j].Hgw-cells[i-1][j].Hgw)/vout;
	if (i < ny-1) if (cells[i+1][j].Hgw < cells[i][j].Hgw) cells[i+1][j].Qstack += cells[i][j].Qstack*cells[i+1][j].Kgw*(cells[i][j].Hgw-cells[i+1][j].Hgw)/vout;
      }	 										       

	/*}*/
 
      cells[i][j].Qw = dx*dy/(spy*dt)*cells[i][j].Qstack;  /*m3 pr. sec*/ 
      cells[i][j].Qstack = 0.0;
      /*cells[i][j].Qw = 0.1*(cells[i][j].Qw - cells[i][j].Qw_old) + cells[i][j].Qw_old;*/

    } 
      else {
	cells[i][j].Qw = 0.0;
	cells[i][j].Qstack = 0.0;
      }

    }


  }/*end calculation of discharge*/ 

}

void glacial_hydrology(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nn = nx*ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double dl = sqrt(dx*dy);

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double A0 = hwprop.A0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double Ls = hwprop.Ls; /*cavity spacing*/ 
  double Lc = hwprop.Lc; /*channel spacing*/
  double alpha = hwprop.alpha;
  double minqw = hwprop.minqw;
  double B = 73.3e6;
  double B0 = B/(pow(spy,1.0/3.0)*rho*g);
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,trans,dhw,Kgw,dSdt,qwpsi;
  double te_s = 0.0;
  double te_c = 0.0;
  double te_ss = 0.0;
  double vout;
  double minice = 10.0;
  double dtw = 0.01; /*artificial time step - y*/
 
  

#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,trans,dhw,te_s,te_c,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,dl,rho,g,rho_w,h0,kh,ds,ss0,B,alpha,kmin,spy,Ls,Lc,Lh,dt,S0,A0,minqw,B0,dtw)
  {


    /*compute steady state efffective pressure in cavities and channels*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) { 
	cells[i][j].qw = cells[i][j].Qw/dl; /*m2 pr. sec*/
	if (cells[i][j].qw < minqw) cells[i][j].qw = minqw;
	cells[i][j].psi = rho*g*(cells[i][j].tslope+1.0e-3); /*kg m-2 s-2*/
	cells[i][j].mw = cells[i][j].qw*cells[i][j].psi/(rho*Lh); /*m s-1*/
 
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
 
	  /*step height*/
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = Ls*sslope;
	  cells[i][j].hs = hs;

	  /*Water pressure*/
	  cells[i][j].Pw = cells[i][j].qw/(cells[i][j].Kgw*cells[i][j].psi+1.0e-16);
	  if (cells[i][j].Pw > 0.9*cells[i][j].tn) cells[i][j].Pw = 0.9*cells[i][j].tn;
	  cells[i][j].te = cells[i][j].tn - cells[i][j].Pw;
  

	  /*effective pressure*/
	  /*cells[i][j].Pw = cells[i][j].tn - B/(rho*g)*pow((Lc*cells[i][j].mw+(cells[i][j].sliding*hs)/spy)/(cells[i][j].Ac+alpha*hs*cells[i][j].S+A0),1.0/3.0);
	  if (cells[i][j].Pw > 0.9*cells[i][j].tn) cells[i][j].Pw = 0.9*cells[i][j].tn;
	  if (cells[i][j].Pw < 0.01*cells[i][j].tn) cells[i][j].Pw = 0.01*cells[i][j].tn;
	  cells[i][j].te = cells[i][j].tn - cells[i][j].Pw;*/


	  /*limit effective pressure*/
	  /*if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	    if (cells[i][j].te < 0.01*cells[i][j].tn) cells[i][j].te = 0.01*cells[i][j].tn;*/

	  /*update cavity size*/
	  cells[i][j].S = (alpha*cells[i][j].S+dtw*cells[i][j].sliding)/(alpha*(1.0+dtw*pow(cells[i][j].te/B0,3.0)));
	  if (cells[i][j].S < S0) cells[i][j].S = S0;
	  if (cells[i][j].S > Ls) cells[i][j].S = Ls;
	  cells[i][j].SLf = cells[i][j].S/Ls;

	  /*update channels*/
	  cells[i][j].Ac = (cells[i][j].Ac+dtw*spy*Lc*cells[i][j].mw)/(1.0+dtw*pow(cells[i][j].te/B0,3.0));
 
	  /*hydrological head*/
          /*cells[i][j].Hgw = cells[i][j].bed + cells[i][j].Pw;*/

	  /*cells[i][j].Kgw = 0.1*kmin;*/

	  /*cells[i][j].Kgw = 0.9*cells[i][j].Kgw + 0.1*kmin/(sqrt(cells[i][j].psi)+0.01)*(pow(alpha*hs*cells[i][j].S,1.2)/Ls);*/

	  cells[i][j].Kgw = 0.9*cells[i][j].Kgw + 0.1*kmin/(sqrt(cells[i][j].psi)+0.01)*(pow(cells[i][j].Ac,1.2)/Lc+pow(alpha*hs*cells[i][j].S,1.2)/Ls);

	}
	else { 
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0;
	  cells[i][j].SLf = 0.0;
	  cells[i][j].hydro = 0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/


void glacial_hydrology_old(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nn = nx*ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double dl = sqrt(dx*dy);

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double A0 = hwprop.A0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double Ls = hwprop.Ls; /*cavity spacing*/ 
  double Lc = hwprop.Lc; /*channel spacing*/
  double alpha = hwprop.alpha;
  double minqw = hwprop.minqw;
  double B = 73.3e6;
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,trans,dhw,Kgw,dSdt,qwpsi;
  double te_s = 0.0;
  double te_c = 0.0;
  double te_ss = 0.0;
  double vout;
  double minice = 10.0;
 
  

#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,trans,dhw,te_s,te_c,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,dl,rho,g,rho_w,h0,kh,ds,ss0,B,alpha,kmin,spy,Ls,Lc,Lh,dt,S0,A0,minqw)
  {


    /*compute steady state efffective pressure in cavities and channels*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) { 
	cells[i][j].qw = cells[i][j].Qw/dl; /*m2 pr. sec*/
	if (cells[i][j].qw < minqw) cells[i][j].qw = minqw;
	cells[i][j].psi = rho*g*(cells[i][j].tslope+1.0e-3); /*kg m-2 s-2*/
	cells[i][j].mw = cells[i][j].qw*cells[i][j].psi/(rho*Lh); /*m s-1*/
 
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
 
	  /*cavity length and pressure*/
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = Ls*sslope;
	  cells[i][j].hs = hs;
  
	  /*cavity size and steady state effective pressure*/
	  temp = cells[i][j].S;
	  cells[i][j].S = pow(Ls*cells[i][j].qw/(kmin*sqrt(cells[i][j].psi)),0.8)/(alpha*hs);
	  if (cells[i][j].S < S0) cells[i][j].S = S0;
	  if (cells[i][j].S > Ls) cells[i][j].S = Ls;
	  cells[i][j].SLf = cells[i][j].S/Ls;
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_s = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S,2.0),1.0/3.0); 

	  /*channel cross section and steady-state effective pressure*/ 
	  temp = cells[i][j].Ac;
	  cells[i][j].Ac = pow(Lc*cells[i][j].qw/(kmin*sqrt(cells[i][j].psi)),0.8);
	  cells[i][j].dAcdt = (cells[i][j].Ac - temp)/dt;
	  temp =  Lc*cells[i][j].mw;
	  te_c = 1.0/(rho*g)*3.0*B*pow(temp/(2.0*(cells[i][j].Ac+A0)),1.0/3.0); 
 
	  /*if cavity drained*/	  
	  if (te_s > te_c) {
	    cells[i][j].te = te_s;
	    cells[i][j].hydro = -1;
	  }
	  else {
	    cells[i][j].te = te_c;
	    cells[i][j].hydro = 1;
	  }

	  /*limit effective pressure*/
	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  if (cells[i][j].te < 0.01*cells[i][j].tn) cells[i][j].te = 0.01*cells[i][j].tn;
 
	  /*hydrological head*/
          cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;

	  /*water pressure*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;

	  /*Transmissivity*/
	  cells[i][j].Kgw = kmin;

	}
	else { 
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0;
	  cells[i][j].SLf = 0.0;
	  cells[i][j].hydro = 0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/



void glacial_hydrology_Iverson_transient(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt,double time)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double L = hwprop.Ls;
  double alpha = hwprop.alpha;
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,dhw,Kgw,dSdt,qwpsi;
  double te_ss = 0.0;


#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,dhw,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,rho,g,rho_w,h0,kh,ds,ss0,minhw,maxhw,tscale,B,alpha,kmin,spy,L,Lh,dt,S0,time)
  {


    /*compute water flux qwx m2/y at time t*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	if (cells[i][j-1].Hgw > cells[i][j].Hgw) {
	  hw = cells[i][j-1].hw; 
	}
	else {
	  hw = cells[i][j].hw;  
	}
	if (hw > 1.0) hw = 1.0; if (hw < 0.0) hw = 0.0;
	hp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i][j-1].Hgw)/dx;
	dhw = -dt*spy*tscale*kmin*pow(L*hw,1.2)*hp[i][j].psi/(sqrt(fabs(hp[i][j].psi)+10.0)*dx*L);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].dhw;
	} 
	else {
	  if (cells[i][j-1].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i][j-1].hw) dhw = cells[i][j-1].hw;
	}
	hp[i][j].dhw = dhw;
      }
    }
    

    /*compute water flux qw m2/yr at time t*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i-1][j].Hgw > cells[i][j].Hgw) { 
	  hw = cells[i-1][j].hw; 
	}
	else {
	  hw = cells[i][j].hw; 
	}
	if (hw < 0.0) hw = 0.0; if (hw > 1.0) hw = 1.0;
	vp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i-1][j].Hgw)/dy;
	dhw = -dt*spy*tscale*kmin*pow(L*hw,1.2)*vp[i][j].psi/(sqrt(fabs(vp[i][j].psi)+10.0)*dy*L);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0; 
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].hw;
	}
	else {
	  if (cells[i-1][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i-1][j].hw) dhw = cells[i-1][j].hw;
	}
	vp[i][j].dhw = dhw;
      }
    }

    /*compute change in water sheet thickness and efffective pressure*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
	  cells[i][j].qw = 0.5*(fabs(hp[i][j].dhw)*dx + fabs(hp[i][j+1].dhw)*dx + fabs(vp[i][j].dhw)*dy + fabs(vp[i+1][j].dhw)*dy)/dt;
	  qwpsi = 0.5*dx*(fabs(hp[i][j].dhw*hp[i][j].psi) + fabs(hp[i][j+1].dhw*hp[i][j+1].psi)) + 0.5*dy*(fabs(vp[i][j].dhw*vp[i][j].psi) + fabs(vp[i+1][j].dhw*vp[i+1][j].psi))/dt;
	  cells[i][j].mw = qwpsi/(tscale*rho*Lh);
	  cells[i][j].dhw = hp[i][j].dhw - hp[i][j+1].dhw + vp[i][j].dhw - vp[i+1][j].dhw + tscale*cells[i][j].water*dt; 
	  cells[i][j].dhwdt = cells[i][j].dhw/dt; 
	  /*if (cells[i][j].dhwdt < -1.0) cells[i][j].dhwdt = -1.0;
	    if (cells[i][j].dhwdt > 1.0) cells[i][j].dhwdt = 1.0;*/
	  cells[i][j].hw += dt*cells[i][j].dhwdt; /*at time t + dt*/
	  if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;
	  if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = L*sslope;
	  cells[i][j].hs = hs;
	  cells[i][j].S = L*cells[i][j].hw/(alpha*hs); 
	  cells[i][j].dSdt = L*cells[i][j].dhwdt/(alpha*hs); 
	  if (cells[i][j].S < 0.0) cells[i][j].S = 0.0;
	  /*if (cells[i][j].S > 2.0) cells[i][j].S = 2.0;*/
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_ss = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  if (te_ss > cells[i][j].tn) te_ss = cells[i][j].tn;
	  temp += L*(cells[i][j].mw-cells[i][j].Mb)/spy;
	  if (time > 0.1) temp -= 0.0*alpha*hs*cells[i][j].dSdt/spy;
	  if (temp > 0.0) {
	    cells[i][j].te_new = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  }
	  else cells[i][j].te_new = 0.01*cells[i][j].tn;
	  if (cells[i][j].te_new < 0.01*cells[i][j].tn) cells[i][j].te_new = 0.01*cells[i][j].tn;
	  if (cells[i][j].te_new > cells[i][j].tn) cells[i][j].te_new = cells[i][j].tn;

	  cells[i][j].te = cells[i][j].te_new;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - te_ss;/*cells[i][j].te;*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  cells[i][j].hw = minhw;
	  cells[i][j].dhwdt = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanhw += cells[i][j].hw;
	meandhwdt += cells[i][j].dhwdt;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanhw /= (double)(nice);
  meandhwdt /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanhw = meanhw;
  (*mesh).meandhwdt = meandhwdt;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/


void glacial_hydrology_Iverson_steady(celltype **cells,meshtype *mesh,iproptype iprop,hwproptype hwprop)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double S0 = hwprop.S0;
  double kh = hwprop.kh;
  double h0 = hwprop.h0;
  double L = hwprop.Ls;
  double alpha = hwprop.alpha;
  double rho_w = 1000.0;
  double k0 = hwprop.kmin; 
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double tb0 = 0.001;
  double Rbed,meante,hw,S,sslope,hs,psi;

#pragma omp parallel shared(cells) private(i,j,Rbed,hw,S,sslope,hs,psi) firstprivate(nx,ny,dx,dy,k0,rho,g,rho_w,h0,kh,minhw,maxhw,tscale,B,spy,alpha,L,tb0)
  {


    /*compute change in water sheet thickness at time t and efffective pressure at time t+dt/2*/
#pragma omp for schedule(static)
    for (i=1;i<(ny-1);i++) {
      for (j=1;j<(nx-1);j++) {
	if (cells[i][j].ice > 10.0) {
	  /*cells[i][j].qw = qfac*cells[i][j].ice*(cells[i][j].deformation+cells[i][j].sliding);*/
	  /*psi = rho*g*(cells[i][j].tslope + tb0); 
	    cells[i][j].hw = pow(cells[i][j].qw*vis_w/(k0*psi),1.0/3.0);*/ 
	  cells[i][j].hw = 0.1*cells[i][j].water;
	  if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;
	  if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;
	  sslope = cells[i][j].slidingslope; if (sslope < 0.0) sslope = 0.0; 
	  hs = kh*sslope + h0; 
	  Rbed = hs/L;
	  S = cells[i][j].hw/(alpha*Rbed); if (S < 1.0e-3) S = 1.0e-3;
	  cells[i][j].te = 1.0/(rho*g)*B*pow((16.0*cells[i][j].sliding/spy*hs/(2.0*pi))/(S*S),1.0/3.0); /*at time t + dt*/ 
	  if (cells[i][j].te < 0.1*cells[i][j].tn) cells[i][j].te = 0.1*cells[i][j].tn;
	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  cells[i][j].hw = minhw;
	  cells[i][j].dhwdt = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.5*cells[i][j].tn;
	}
      }
    } 

#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      cells[i][0].te = cells[i][1].te;
      cells[i][nx-1].te = cells[i][nx-2].te;
    }

#pragma omp for schedule(static)
    for (j=0;j<nx;j++) {
      cells[0][j].te = cells[1][j].te;
      cells[ny-1][j].te = cells[ny-2][j].te;
    }

  }/*pragma*/

  /*compute mean effective pressure*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meante += cells[i][j].te;
    }
  }
  meante /= (double)((*mesh).nc);
  (*mesh).meante = meante;

}/*void*/


void glacial_hydrology_transient(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt,double time)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double A0 = hwprop.A0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double Ls = hwprop.Ls; /*cavity spacing*/ 
  double Lc = hwprop.Lc; /*channel spacing*/
  double alpha = hwprop.alpha;
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,dhw,Kgw,dSdt,qwpsi;
  double te_s = 0.0;
  double te_c = 0.0;
  double te_ss = 0.0;

#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,dhw,te_s,te_c,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,rho,g,rho_w,h0,kh,ds,ss0,minhw,maxhw,tscale,B,alpha,kmin,spy,Ls,Lc,Lh,dt,S0,A0,time)
  {


    /*compute water flux qwx m2/y at time t*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	if (cells[i][j-1].Hgw > cells[i][j].Hgw) {
	  hw = cells[i][j-1].hw; 
	}
	else {
	  hw = cells[i][j].hw;  
	}
	if (hw > 1.0) hw = 1.0; if (hw < 0.0) hw = 0.0;
	hp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i][j-1].Hgw)/dx;
	dhw = -dt*spy*tscale*kmin*pow(Ls*hw,1.2)*hp[i][j].psi/(sqrt(fabs(hp[i][j].psi)+10.0)*dx*Ls);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].dhw;
	} 
	else {
	  if (cells[i][j-1].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i][j-1].hw) dhw = cells[i][j-1].hw;
	}
	hp[i][j].dhw = dhw;
      }
    }
    

    /*compute water flux qw m2/yr at time t*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i-1][j].Hgw > cells[i][j].Hgw) { 
	  hw = cells[i-1][j].hw; 
	}
	else {
	  hw = cells[i][j].hw; 
	}
	if (hw < 0.0) hw = 0.0; if (hw > 1.0) hw = 1.0;
	vp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i-1][j].Hgw)/dy;
	dhw = -dt*spy*tscale*kmin*pow(Ls*hw,1.2)*vp[i][j].psi/(sqrt(fabs(vp[i][j].psi)+10.0)*dy*Ls);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0; 
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].hw;
	}
	else {
	  if (cells[i-1][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i-1][j].hw) dhw = cells[i-1][j].hw;
	}
	vp[i][j].dhw = dhw;
      }
    }

    /*compute change in water sheet thickness and efffective pressure*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].qw = 0.5*(fabs(hp[i][j].dhw + hp[i][j+1].dhw)*dx + fabs(vp[i][j].dhw + vp[i+1][j].dhw)*dy)/dt;
	qwpsi = 0.5*dx*(fabs(hp[i][j].dhw*hp[i][j].psi) + fabs(hp[i][j+1].dhw*hp[i][j+1].psi)) + 0.5*dy*(fabs(vp[i][j].dhw*vp[i][j].psi) + fabs(vp[i+1][j].dhw*vp[i+1][j].psi))/dt;
	cells[i][j].mw = qwpsi/(tscale*rho*Lh);
	cells[i][j].dhw = hp[i][j].dhw - hp[i][j+1].dhw + vp[i][j].dhw - vp[i+1][j].dhw + tscale*cells[i][j].water*dt; 
	cells[i][j].dhwdt = cells[i][j].dhw/dt; 
	/*if (cells[i][j].dhwdt < -1.0) cells[i][j].dhwdt = -1.0;
	  if (cells[i][j].dhwdt > 1.0) cells[i][j].dhwdt = 1.0;*/
	cells[i][j].hw += dt*cells[i][j].dhwdt; /*at time t + dt*/
	/*if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;*/
	/*if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;*/
 
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
 
	  /*cavity length and pressure*/
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = Ls*sslope;
	  cells[i][j].hs = hs;
	  cells[i][j].S = Ls*cells[i][j].hw/(alpha*hs); 
	  cells[i][j].dSdt = Ls*cells[i][j].dhwdt/(alpha*hs); 
	  if (cells[i][j].S < 0.0) cells[i][j].S = 0.0;
	  /*if (cells[i][j].S > 2.0) cells[i][j].S = 2.0;*/
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_ss = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  /*temp -= alpha*hs*cells[i][j].dSdt/spy;*/
	  te_s = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 

	  /*channel cross section and pressure*/
	  cells[i][j].Ac = Lc*cells[i][j].hw;
	  cells[i][j].dAcdt = Lc*cells[i][j].dhwdt;
	  temp =  Lc*(cells[i][j].mw-cells[i][j].Mb)/spy;
	  /*temp -= cells[i][j].dAcdt/spy;*/
	  te_c = 1.0/(rho*g)*3.0*B*pow(temp/(2.0*(cells[i][j].Ac+A0)),1.0/3.0); 
 
	  /*if cavity drained*/	  
	  if (te_s > te_c) {
	    cells[i][j].te = te_s;
	    cells[i][j].hydro = -1;
	  }
	  else {
	    cells[i][j].te = te_c;
	    cells[i][j].hydro = 1;
	  }

	  if (te_ss > cells[i][j].tn) te_ss = cells[i][j].tn;
	  if (te_ss < 0.01*cells[i][j].tn) te_ss = 0.01*cells[i][j].tn;


	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  if (cells[i][j].te < 0.01*cells[i][j].tn) cells[i][j].te = 0.01*cells[i][j].tn;
 
	  cells[i][j].te_s = te_ss;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - te_ss;
          /*cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  /*cells[i][j].hw = minhw;*/
	  /*cells[i][j].dhwdt = 0.0;*/
	  cells[i][j].Hgw = cells[i][j].bed + 1.0e5*cells[i][j].hw;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0; 
	  cells[i][j].hydro = 0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanhw += cells[i][j].hw;
	meandhwdt += cells[i][j].dhwdt;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanhw /= (double)(nice);
  meandhwdt /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanhw = meanhw;
  (*mesh).meandhwdt = meandhwdt;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/
