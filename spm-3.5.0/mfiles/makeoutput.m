function makeoutput(fnr)

SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
L = SPM.mesh.L;
H = SPM.mesh.H;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;

loaddata; 

%make mesh
[X,Y] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
  
%slope perpendicular erosion
erosion = weathering+hillslope+periglacial; 
 
mdata.fnr = fnr;
mdata.X = X;
mdata.Y = Y;
mdata.bed = bed;
mdata.ice = ice;
mdata.sliding = sliding;
mdata.deformation = deformation;
mdata.erosion = erosion;
mdata.sediment = sediment;
mdata.ssedi = ssedi;
mdata.massb = massb;

%load particles
zp = 0;
loadparticles;
pdata.xp = xp;
pdata.yp = yp;
pdata.zp = zp;
pdata.bf = bf;
pdata.sedi = sedi;
pdata.N10 = N10;
pdata.bx = bx;
pdata.by = by;
pdata.birthday = birthday;
pdata.age = age;
pdata.dl = dl;
pdata.vx = vx;
pdata.vy = vy;
pdata.vz = vz;

save(['./output/output',num2str(fnr),'.mat'],'mdata','pdata');