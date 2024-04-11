function prepare_flowlines(fnr)

%load input
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
L = SPM.mesh.L;
H = SPM.mesh.H;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;

%load output
loaddata;

%save data needed for flowline calculations
save data_for_flowlines.mat SPM bed ice sediment icemargin vx vy ...
    massb bmelt