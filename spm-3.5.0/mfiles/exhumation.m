function exhumation(fnrs)

close all;

SPM = SPMload;
L = SPM.mesh.L;
H = SPM.mesh.H;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

Nt = length(fnrs);
Nc = ny*nx;

ex = zeros(Nc,Nt);
burial = zeros(Nc,Nt);
icehist = zeros(Nc,Nt);
topohist = zeros(Nc,Nt);

h = waitbar(0,'Reading files...');
%accumulated erosion
for kk=1:Nt,
    
    fnr = fnrs(kk);
    loaddata;
    
    icehist(:,kk) = ice(:);
    ex(:,kk) = abrasion(:);
    topohist(:,kk) = bed(:);
    
    waitbar(kk/Nt,h);
end; 
close(h);

%burial depth
for i=1:Nt,
    burial(:,i) = ex(:,Nt) - ex(:,i);
end;

save burial.mat ex burial icehist topohist;

surf(Xc,Yc,reshape(burial(:,Nt-1),[ny,nx]));