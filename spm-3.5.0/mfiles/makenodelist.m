function makenodelist()

SPM = SPMload;

nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;


%mouse coordinate input
[px,py]=ginput;

for i=1:length(px),
    
    nnx(i) = floor(px(i)/dx);
    nny(i) = floor(py(i)/dy);
    
    if nnx(i) < 1, nnx(i) = 1; end;
    if nny(i) < 1, nny(i) = 1; end;
    if nnx(i) > nx, nnx(i) = nx; end;
    if nny(i) > ny, nny(i) = ny; end;
    
end;
nnx = nnx(:);
nny = nny(:);

save ./input/nodelist.mat nnx nny;
