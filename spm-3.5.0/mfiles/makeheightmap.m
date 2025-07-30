function makeheightmap(fnr)

close all;

%load SPM
SPM = SPMload;
L = SPM.mesh.L
H = SPM.mesh.H
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
dx = L/nx;
dy = H/ny;

%bitmap matrices
xx = linspace(0,L,nx);
yy = linspace(0,H,ny);
[Xm,Ym] = meshgrid(xx,yy);

%compute dimensions
%dd = 10;
dd = min(dx,dy);
hnx = ceil(L/dd)
hny = ceil(H/dd)

%bitmap matrices
xx = linspace(0,L,hnx);
yy = linspace(0,H,hny);
[Xh,Yh] = meshgrid(xx,yy);


di = dd/5;
inx = ceil(L/di)
iny = ceil(H/di)

%bitmap matrices
xx = linspace(0,L,inx);
yy = linspace(0,H,iny);
[Xi,Yi] = meshgrid(xx,yy);

if fnr == -1,
    
  bed = SPM.data.bed;
  ice = SPM.data.ice;
  
else,
  
  loaddata;
        
end;

%smooth ice thickness
%for i=1:5,
%  ice = smooth(ice);
%end;
  
%interpolate topography
Z = interp2(Xm,Ym,bed+sediment+ice,Xh,Yh);

%ice mask
ICE = interp2(Xm,Ym,ice,Xh,Yh);
Zim = zeros(size(ICE));
I = find(ICE > 5);
Zim(I) = 1;

%corner grid
xx = linspace(-dd/2,L+dd/2,hnx+1);
yy = linspace(-dd/2,H+dd/2,hny+1);
[Xc,Yc] = meshgrid(xx,yy);

Cim = zeros(hny+1,hnx+1);
Cim(2:end-1,2:end-1) = .25*(Zim(1:end-1,1:end-1)+Zim(2:end,1:end-1)+Zim(2:end,2:end)+Zim(1:end-1,2:end));

I = find(Cim <= 0.6); Cim(I) = 0;
I = find(Cim > 0.6); Cim(I) = 1;

Zimi = interp2(Xc,Yc,Cim,Xi,Yi);


%for i=1:5,
%  Zim = smooth(Zim);
%end;

%I = find(Zim < 0.95);
%Zim(I) = 0;


%sediment mask
SEDI = interp2(Xm,Ym,sediment,Xh,Yh);
SEDIm = zeros(size(Z));
I = find(SEDI > 20);
SEDIm(I) = 1;

%set colormap
val = linspace(0,1,2^8);
map = val'*ones(1,3);
Nc = length(map(:,1));
colormap(map);

%scale matrices
Zimi = Nc*Zimi;
SEDIm = Nc*SEDIm;

%plots
if 1,
subplot(2,2,1);
imagesc(Z); set(gca,'ydir','normal');
subplot(2,2,4);
imagesc(Zimi); set(gca,'ydir','normal');
subplot(2,2,3);
imagesc(SEDIm); set(gca,'ydir','normal');
end



maketer(Z,dd,fnr);
%maketer(ICE,dd,fnr);
imwrite(flipud(Zimi),map,['./hmaps/icemask',num2str(fnr),'.bmp'],'bmp');
imwrite(flipud(SEDIm),map,['./hmaps/sedimask',num2str(fnr),'.bmp'],'bmp');

