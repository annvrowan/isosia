
void update_particles(celltype **cells,hptype **hp,vptype **vp,meshtype mesh,parproptype parprop,hproptype hprop,particletype *pp,long npmax,long *np,long *npa,long *npia,long *pactive,long *pinactive,long **pindex,double time,double dt)
{

  int ni,nj,kill,ranp;
  long i,j,ip,rip,cnr;
  double xxn,yyn,zzn;
  double vxx,vyy,vzz; 
  double vxx_d,vyy_d,vxx_b,vyy_b; 
  double prate,ran,Kdiff,fac,tfac,ddl;

  double dx = mesh.dx;
  double dy = mesh.dy;
  int nx = mesh.nx;
  int ny = mesh.ny;
  int dodeposit = mesh.dodeposit;
  int dohillslope = mesh.dohillslope;

  /*Initiate*/
  long ia = -1; /*index for active particles*/ 
  double minice = parprop.minice;
  double minsedi = parprop.minsedi; 
  double maxsedi = parprop.maxsedi;
  int maxp = parprop.maxp; 
  int maxpm = parprop.maxpm;
  int minpm = parprop.minpm; 
  double minsedim = parprop.minsedim; 

  /*Cosmo params*/ 
  double erate = 1.0e-1; /*steady state erosion rate cm/yr*/
  double Lspal10 = 150.0; /*Attennuation length spallation 10Be g/cm2 - input instead*/ 
  double lambda10 = 5.0301e-7; /*Decay constant 10Be 1/yr*/
  double rhoi = 0.9; /*density of ice g/cm3*/ 
  double rhor = 2.65;/*density of rock g/cm3*/
  double speed;
  double maxspeed = .25*(dx+dy)/dt; 
  double Ks = hprop.Ks;
  double sc = hprop.sc; 

  double L0 = 3000.0;
  double f0 = 1.0;

  /*loop cells to form particles from sediment*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
 
      /*criteria for particle formation*/
      /*if (((cells[i][j].ice > minice)&&(cells[i][j].sedi > minsedi)&&(cells[i][j].np < maxp))||((cells[i][j].margin > 0)&&(cells[i][j].np < minpm)&&(minsedim > 0.0))) {*/
      /*if (((cells[i][j].ice > minice)&&(cells[i][j].sedi > minsedi)&&(cells[i][j].np < maxp))) {*/
      if (((cells[i][j].sedi > minsedi)&&(cells[i][j].np < maxp))) {

	/*get particle from inactive stack*/
	if ((*npia) >= 0) {
	  ip = pinactive[(*npia)];
	  (*npia) -= 1;
	}
	/*else increase number of particles*/
	else if ((*np) < npmax-1) {
	  (*np) += 1;
	  ip = (*np);
	}
	/*or do not form*/
	else ip = -1;

	/*form particle*/
	if (ip >= 0) {
	  pp[ip].active = 1;
	  ran = (double)rand()/(double)RAND_MAX;
	  pp[ip].x = cells[i][j].x+dx*(ran-0.5);
	  ran = (double)rand()/(double)RAND_MAX;
	  pp[ip].y = cells[i][j].y+dy*(ran-0.5);
	  pp[ip].birthday = time;
	  pp[ip].age = 0.0;
	  pp[ip].bx = pp[ip].x;
	  pp[ip].by = pp[ip].y;
	  pp[ip].dl = 0.0; 
	  pp[ip].afac = f0;

	  /*erosion rate in cell cm/y*/
	  erate = 100.0*(cells[i][j].hillslope_erate+cells[i][j].weathering_rate+cells[i][j].periglacial_erate);

	  /*pickup sediment*/
	  if (cells[i][j].sedi < maxsedi) {
	    pp[ip].sedi = cells[i][j].sedi;
	    cells[i][j].sedi = 0.0;  
	  }
	  else {
	    pp[ip].sedi = maxsedi; 
	    cells[i][j].sedi -= maxsedi; 
	  }


	  /*initialize burial depth and N10*/
	  if (cells[i][j].ice > minice) {  
	    pp[ip].bf = 1.0;  /*sediment is formed under ice*/  
	    pp[ip].N10 = cells[i][j].CNprod/(lambda10+rhor*erate/cells[i][j].atten); /*!!!!!!!*/
	  }
	  else {
	    pp[ip].bf = 0.0; /*noburial*/ 
	    pp[ip].N10 = cells[i][j].CNprod/(lambda10+rhor*erate/cells[i][j].atten);	  
	  }
	  pp[ip].erate = erate; /*save erate on particle*/
	  (*npa) += 1;
	  pactive[(*npa)] = ip;

	}
	  
      }

      /*reset surface and bed sediment arrarys*/
      cells[i][j].ssedi = 0.0;
      cells[i][j].bsedi = 0.0;
      cells[i][j].msedi = 0.0;  
      cells[i][j].afac = 0.0;
      cells[i][j].sN10 = 0.0;
      cells[i][j].np = 0;
    }
  }
  
  
  /*Loop active particles to update positions and kill*/
  for (i=0;i<(*npa);i++) {

    ip = pactive[i]; /*particle number*/  
    
    /*find cell number of particle*/
    ni = (int)floor(pp[ip].y/dy);
    nj = (int)floor(pp[ip].x/dx);

    if (ni < 0) ni = 0;
    else if (ni > ny-1) ni = ny-1;
    if (nj < 0) nj = 0;
    else if (nj > nx-1) nj = nx-1;
    
    /*default: do not kill*/
    kill = 0;

    /*if particle qualifies - particles along mesh boundary dies*/    
    if ((ni <= 0)||(ni >= ny-1)||(nj <= 0)||(nj >= nx-1)||(pp[ip].active == 0)) {
      kill = 1;
    }
    else if (cells[ni][nj].pit == 1) kill = 1;
    else {
      if (cells[ni][nj].np > maxp) { /*merge*/
	kill = 2;
      }
    }


    if (kill == 0) {

      /*count cell*/
      cells[ni][nj].np += 1;

      /*global cell number*/
      cnr = ni + nj*ny;

      /*registre particle*/
      pindex[cnr][cells[ni][nj].np-1] = ip;

      /*if ice covered*/
      if (cells[ni][nj].ice > minice) {

	/*compute velocity using bilinear interpolation*/
	yyn = (pp[ip].y-(double)ni*dy)/dy;
	xxn = (pp[ip].x-(double)nj*dx)/dx;
	vxx_d = (1.0-xxn)*(1.0-yyn)*hp[ni][nj].vx_d+(1.0-xxn)*yyn*hp[ni+1][nj].vx_d+xxn*(1.0-yyn)*hp[ni][nj+1].vx_d+xxn*yyn*hp[ni+1][nj+1].vx_d;
	vyy_d = (1.0-xxn)*(1.0-yyn)*vp[ni][nj].vy_d+(1.0-xxn)*yyn*vp[ni+1][nj].vy_d+xxn*(1.0-yyn)*vp[ni][nj+1].vy_d+xxn*yyn*vp[ni+1][nj+1].vy_d;
	vxx_b = (1.0-xxn)*(1.0-yyn)*hp[ni][nj].vx_b+(1.0-xxn)*yyn*hp[ni+1][nj].vx_b+xxn*(1.0-yyn)*hp[ni][nj+1].vx_b+xxn*yyn*hp[ni+1][nj+1].vx_b;
	vyy_b = (1.0-xxn)*(1.0-yyn)*vp[ni][nj].vy_b+(1.0-xxn)*yyn*vp[ni+1][nj].vy_b+xxn*(1.0-yyn)*vp[ni][nj+1].vy_b+xxn*yyn*vp[ni+1][nj+1].vy_b;

	vxx = 1.2*(1.0-pow(pp[ip].bf,4.0))*vxx_d + vxx_b;
	vyy = 1.2*(1.0-pow(pp[ip].bf,4.0))*vyy_d + vyy_b;
	vzz = (1.0-pp[ip].bf)*cells[ni][nj].Ms-pp[ip].bf*cells[ni][nj].Mb; /*positive downward - note that Mb is negative*/
  
	/*update vertical position*/
	pp[ip].bf += vzz*dt/cells[ni][nj].ice;
	if (pp[ip].bf < 0.0) pp[ip].bf = 0.0;
	else if (pp[ip].bf > 1.0) pp[ip].bf = 1.0;
  
      }
      else if (dohillslope > 0) {
 

	fac = 1.0 - pow(cells[ni][nj].bslope/sc,2.0); if (fac < 0.001) fac = 0.001;
	Kdiff = Ks/fac;

	vxx = -Kdiff*cells[ni][nj].dbdx;
	vyy = -Kdiff*cells[ni][nj].dbdy; 
	vzz = 0.0; 
	pp[ip].bf = 0.0; /*this causes perticles that move onto ice, to start at the ice surface*/

	speed = sqrt(vxx*vxx+vyy*vyy);
	if (speed >= maxspeed) {
	  vxx *= maxspeed/speed;
	  vyy *= maxspeed/speed;  
	}

      }

      else {

	vxx = 0.0;
	vyy = 0.0;
	vzz = 0.0;

      }


      /*store velocity*/
      pp[ip].vx = vxx;
      pp[ip].vy = vyy;
      pp[ip].vz = vzz;

      /*update horizontal positions*/
      pp[ip].x += vxx*dt;
      pp[ip].y += vyy*dt;

      /*vertical position*/
      pp[ip].z = cells[ni][nj].bed + (1.0-pp[ip].bf)*cells[ni][nj].ice; 
 
      /*update CNs*/  
      prate = cells[ni][nj].CNprod*exp(-rhoi*pp[ip].bf*cells[ni][nj].ice*100.0/cells[ni][nj].atten);
      pp[ip].N10 = pp[ip].N10*exp(-dt*lambda10)+prate*(1.0-exp(-dt*lambda10))/lambda10;

      /*update other things*/
      pp[ip].age += dt; 
      ddl = dt*sqrt(vxx*vxx+vyy*vyy);
      pp[ip].dl += ddl;
      tfac = (2.0*L0-ddl)/(2.0*L0+ddl); if (tfac < 0.0) tfac = 0.0;
      pp[ip].afac *= tfac; 

      /*add to surface and bed sediment arrays*/ 
      if (pp[ip].bf < 0.05) cells[ni][nj].ssedi += pp[ip].sedi; 
      else if (pp[ip].bf > 0.95) {
	cells[ni][nj].bsedi += pp[ip].sedi; 
	cells[ni][nj].afac += pp[ip].sedi*pp[ip].afac; 
      }
      else  cells[ni][nj].msedi += pp[ip].sedi; 
 
      /*add to average N10 of surface sediment*/ 
      if (pp[ip].bf < 0.05) cells[ni][nj].sN10 += pp[ip].sedi*pp[ip].N10; 


    }

    /*merge particle with another in the cell*/
    else if (kill == 2) {

      /*choose random particle*/
      ran = (double)rand()/(double)RAND_MAX;
      ranp = (int)floor((double)ran*(double)cells[ni][nj].np);
      if (ranp < 0) ranp = 0;
      if (ranp > cells[ni][nj].np-1) ranp = cells[ni][nj].np-1;

      /*global cell number*/
      cnr = ni + nj*ny;
      
      /*random partile number*/
      rip = pindex[cnr][ranp];

      pp[rip].N10 = (pp[rip].sedi*pp[rip].N10+pp[ip].sedi*pp[ip].N10)/(pp[rip].sedi+pp[ip].sedi);
      pp[rip].sedi += pp[ip].sedi;
      pp[ip].sedi = 0.0; 
      pp[ip].N10 = 0.0;
      pp[ip].active = 0;
      (*npia) += 1;
      pinactive[(*npia)] = ip;


    }



    /*kill active particle*/
    else {

      /*deposit sediment at ice margin*/      
      /*if (cells[ni][nj].margin >= 0) cells[ni][nj].sedi += pp[ip].sedi;*/ 
      if (dodeposit > 0) cells[ni][nj].sedi += pp[ip].sedi;
      pp[ip].sedi = 0.0; 
      pp[ip].N10 = 0.0;
      pp[ip].active = 0;
      (*npia) += 1;
      pinactive[(*npia)] = ip;

      /*write to file*/

    }
  }
  
 

  /*condence active stack*/
  for (i=0;i<=(*npa);i++) {
    
    ip = pactive[i]; /*particle number*/  

    if (pp[ip].active > 0) {
      ia += 1;
      pactive[ia] = ip;

    }

  }
  (*npa) = ia; /*condence active stack*/
 

}/*update_particles*/
