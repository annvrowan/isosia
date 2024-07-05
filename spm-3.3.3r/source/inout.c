
meshtype read_mesh(char *path,char file[200],FILE *fp)
{

  meshtype mesh;

  /*read from file mesh.input*/
  sprintf(file,"%s/input/mesh.input",path); /*filename*/
  if ((fp = fopen(file,"rt")) == NULL) {
    printf("Cannot open file : %s/input/mesh.input\n",path);
    exit(1);
  }
  fscanf(fp,"%lf\n",&mesh.L); 
  fscanf(fp,"%lf\n",&mesh.H); 
  fscanf(fp,"%d\n",&mesh.nx); 
  fscanf(fp,"%d\n",&mesh.ny); 
  fscanf(fp,"%lf\n",&mesh.dx); 
  fscanf(fp,"%lf\n",&mesh.dy); 
  fscanf(fp,"%lf\n",&mesh.hmin); 
  fscanf(fp,"%lf\n",&mesh.ct); 
  fscanf(fp,"%lf\n",&mesh.maxtime); 
  fscanf(fp,"%lf\n",&mesh.filetime); 
  fscanf(fp,"%lf\n",&mesh.maxdt); 
  fscanf(fp,"%lf\n",&mesh.gravity); 
  fscanf(fp,"%lf\n",&mesh.maxb); 
  fscanf(fp,"%lf\n",&mesh.maxs); 
  fscanf(fp,"%d\n",&mesh.nmoni); 
  fscanf(fp,"%d\n",&mesh.periodic); 
  fscanf(fp,"%d\n",&mesh.slidingmode); 
  fscanf(fp,"%d\n",&mesh.coldbased); 
  fscanf(fp,"%d\n",&mesh.hydromode); 
  fscanf(fp,"%d\n",&mesh.dofluvial); 
  fscanf(fp,"%d\n",&mesh.doice); 
  fscanf(fp,"%d\n",&mesh.dosliding); 
  fscanf(fp,"%d\n",&mesh.doperiglacial); 
  fscanf(fp,"%d\n",&mesh.doglacialhydrology); 
  fscanf(fp,"%d\n",&mesh.doglacialerosion); 
  fscanf(fp,"%d\n",&mesh.dohillslope); 
  fscanf(fp,"%d\n",&mesh.dohillslopeerosion); 
  fscanf(fp,"%d\n",&mesh.dolandslide); 
  fscanf(fp,"%d\n",&mesh.doavalance); 
  fscanf(fp,"%d\n",&mesh.doisostasy); 
  fscanf(fp,"%d\n",&mesh.docelldata); 
  fscanf(fp,"%d\n",&mesh.dosediment); 
  fscanf(fp,"%d\n",&mesh.doglacialsedi); 
  fscanf(fp,"%d\n",&mesh.dodebrisablation); 
  fclose(fp);

  /*number of cells*/
  mesh.nc = mesh.nx*mesh.ny;

  return(mesh);

}

iproptype read_iprop(char *path,char file[200],FILE *fp)
{

  iproptype iprop;

  /*read from file mesh.input*/
  sprintf(file,"%s/input/iprop.input",path); /*filename*/
  if ((fp = fopen(file,"rt")) == NULL) {
    printf("Cannot open file : %s/input/iprop.input\n",path);
    exit(1);
  }
  fscanf(fp,"%lf\n",&iprop.gamma); 
  fscanf(fp,"%lf\n",&iprop.gamma0); 
  fscanf(fp,"%lf\n",&iprop.Cs); 
  fscanf(fp,"%lf\n",&iprop.latentheat); 
  fscanf(fp,"%lf\n",&iprop.ki); 
  fscanf(fp,"%lf\n",&iprop.rho); 
  fscanf(fp,"%lf\n",&iprop.cp); 
  fscanf(fp,"%lf\n",&iprop.ifac); 
  fscanf(fp,"%lf\n",&iprop.sfac); 
  fscanf(fp,"%lf\n",&iprop.vbfac); 
  fscanf(fp,"%d\n",&iprop.maxitt_v); 
  fscanf(fp,"%d\n",&iprop.maxitt_s); 
  fscanf(fp,"%lf\n",&iprop.C); 
  fscanf(fp,"%lf\n",&iprop.L0); 
  fscanf(fp,"%lf\n",&iprop.minbeta); 
  fscanf(fp,"%lf\n",&iprop.maxsliding); 
  fscanf(fp,"%lf\n",&iprop.maxdeformation); 
  fscanf(fp,"%lf\n",&iprop.mf); 
  fscanf(fp,"%lf\n",&iprop.qfac); 
  fscanf(fp,"%lf\n",&iprop.Ka); 
  fscanf(fp,"%lf\n",&iprop.ap); 
  fclose(fp);

  return(iprop);

}

mproptype read_mprop(char *path,char file[200],FILE *fp)
{

  int i,j;
  mproptype mprop;


  /*read from file mprop.input*/
  sprintf(file,"%s/input/mprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/mprop.input\n",path);
    exit(1);
  }
  fread(&mprop.mtype,sizeof(int),1,fp);
  fread(&mprop.avaslope,sizeof(double),1,fp);
  fread(&mprop.avacurv,sizeof(double),1,fp);
  fread(&mprop.maxslope,sizeof(double),1,fp);
  fread(&mprop.lrate,sizeof(double),1,fp);
  fread(&mprop.Tsl,sizeof(double),1,fp);
  fread(&mprop.maxacc,sizeof(double),1,fp);
  fread(&mprop.maxabla,sizeof(double),1,fp);
  fread(&mprop.accgrad,sizeof(double),1,fp);
  fread(&mprop.ablgrad,sizeof(double),1,fp);
  fread(&mprop.qb,sizeof(double),1,fp);
  fread(&mprop.sedifrac,sizeof(double),1,fp);
  fread(&mprop.Ldebris,sizeof(double),1,fp);
  fread(&mprop.nTemp,sizeof(int),1,fp);
  mprop.Temp = malloc(2*sizeof(double*)); for (i=0;i<2;i++) mprop.Temp[i] = malloc(mprop.nTemp*sizeof(double));
  for (j=0;j<mprop.nTemp;j++) for (i=0;i<2;i++) fread(&mprop.Temp[i][j],sizeof(double),1,fp);
  fread(&mprop.nMrate,sizeof(int),1,fp);
  mprop.Mrate = malloc(2*sizeof(double*)); for (i=0;i<2;i++) mprop.Mrate[i] = malloc(mprop.nMrate*sizeof(double));
  for (j=0;j<mprop.nMrate;j++) for (i=0;i<2;i++) fread(&mprop.Mrate[i][j],sizeof(double),1,fp);
  fclose(fp);

  return(mprop);

}

hwproptype read_hwprop(char *path,char file[200],FILE *fp)
{

  hwproptype hwprop;


  /*read from file hwprop.input*/
  sprintf(file,"%s/input/hwprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/hwprop.input\n",path);
    exit(1);
  }
  fread(&hwprop.Kgw,sizeof(double),1,fp);
  fread(&hwprop.a2w,sizeof(double),1,fp);
  fread(&hwprop.po,sizeof(double),1,fp);
  fread(&hwprop.tscale,sizeof(double),1,fp);
  fclose(fp);

  return(hwprop);

}

hproptype read_hprop(char *path,char file[200],FILE *fp)
{

  hproptype hprop;


  /*read from file hprop.input*/
  sprintf(file,"%s/input/hprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/hprop.input\n",path);
    exit(1);
  }
  fread(&hprop.Ks,sizeof(double),1,fp);
  fread(&hprop.sc,sizeof(double),1,fp);
  fread(&hprop.Ke,sizeof(double),1,fp);
  fread(&hprop.Ls,sizeof(double),1,fp);
  fread(&hprop.gamma,sizeof(double),1,fp);
  fread(&hprop.Nc,sizeof(int),1,fp);
  fclose(fp);

  return(hprop);

}

fproptype read_fprop(char *path,char file[200],FILE *fp)
{

  fproptype fprop;


  /*read from file hprop.input*/
  sprintf(file,"%s/input/fprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/fprop.input\n",path);
    exit(1);
  }
  fread(&fprop.pr,sizeof(double),1,fp);
  fread(&fprop.rho_s,sizeof(double),1,fp);
  fread(&fprop.Dg,sizeof(double),1,fp);
  fread(&fprop.kw,sizeof(double),1,fp);
  fread(&fprop.tau_c,sizeof(double),1,fp);
  fread(&fprop.Kt,sizeof(double),1,fp);
  fread(&fprop.Ke,sizeof(double),1,fp);

  fclose(fp);

  return(fprop);

}

pgproptype read_pgprop(char *path,char file[200],FILE *fp)
{

  int i,j;

  pgproptype pgprop;

  /*read periglacial properties*/
  sprintf(file,"%s/input/periglacial.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/periglacial.input\n",path);
    exit(1);
  }
  fread(&pgprop.nHs,sizeof(int),1,fp);
  fread(&pgprop.nT0,sizeof(int),1,fp);

  /*allocate*/
  pgprop.Hsv = malloc(pgprop.nHs*sizeof(double)); 
  pgprop.T0v = malloc(pgprop.nT0*sizeof(double)); 
  pgprop.Ci = malloc(pgprop.nHs*sizeof(double*)); for (i=0;i<pgprop.nHs;i++) pgprop.Ci[i] = malloc(pgprop.nT0*sizeof(double));
  pgprop.Tr = malloc(pgprop.nHs*sizeof(double*)); for (i=0;i<pgprop.nHs;i++) pgprop.Tr[i] = malloc(pgprop.nT0*sizeof(double));

  for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Hsv[i],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) fread(&pgprop.T0v[j],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Ci[i][j],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Tr[i][j],sizeof(double),1,fp);
  fread(&pgprop.rho_b,sizeof(double),1,fp);
  fread(&pgprop.rho_s,sizeof(double),1,fp);
  fread(&pgprop.Ke,sizeof(double),1,fp);
  fread(&pgprop.Kt,sizeof(double),1,fp);
  fread(&pgprop.maxsedi,sizeof(double),1,fp);
  fread(&pgprop.maxice,sizeof(double),1,fp);
  fread(&pgprop.minslope,sizeof(double),1,fp);
  fread(&pgprop.minsedi,sizeof(double),1,fp);
  fclose(fp);

  return(pgprop);

}


void write_output(celltype **cells,hptype **hp,vptype **vp,meshtype mesh,char *path,char file[200],FILE *fp,long step)
{

  int i,j,k;
  double vx,vy;
  int nx = mesh.nx;
  int ny = mesh.ny;

  /*create output file*/
  sprintf(file,"%s/output/output%ld.dat",path,step);
  if ((fp = fopen(file,"wb")) == NULL) {
    printf("could not open file for writing\n");
    exit(1);
  }
  else {
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].x,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].y,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].bed,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].ice,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].bslope,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].tslope,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].te2,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].periglacial_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sedi,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].tn,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].te,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].ts,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sliding,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].deformation,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Pw,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Mb,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].quarrying,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].accrate,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ms,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].periglacial_erate,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].fluvial_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].landslide_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].abrasion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].isostasy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].lee,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ts,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].hillslope_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ta,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Tb,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sfac,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].margin,sizeof(int),1,fp);
    for (k=0;k<20;k++) for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Vs[k],sizeof(double),1,fp);
    fclose(fp);
  }
    
}/*write_output*/
