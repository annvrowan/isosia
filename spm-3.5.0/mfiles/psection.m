function zp=psection(fnr,jpg)


%plot average particle properties in flowline section

close all
set(gcf,'units','normalized','position',[.1,.2,.6,.6]);

pcol1 = [.8,.4,.4];
pcol2 = [.5,.5,.5];

%load moraine positions
load mmoraine.mat;

%load particles data
loadparticles;

%**** form moraine polygon **********
x0 = lx(1);
y0 = ly(1);
dw = 100; %width of moraine

pp1x = lx;
pp1y = ly;
pp2x = lx;
pp2y = ly;

for i=2:length(lx),
    dx = lx(i)-lx(i-1);
    dy = ly(i)-ly(i-1);
    dd = sqrt(dx*dx+dy*dy);
    nx = -dy/dd;
    ny = dx/dd;
    pp1x(i-1) = lx(i-1) + dw*nx;
    pp1y(i-1) = ly(i-1) + dw*ny;
    pp2x(i-1) = lx(i-1) - dw*nx;
    pp2y(i-1) = ly(i-1) - dw*ny;
end;
pp1x(end) = lx(end) + dw*nx;
pp1y(end) = ly(end) + dw*ny;
pp2x(end) = lx(end) - dw*nx;
pp2y(end) = ly(end) - dw*ny;

ppx = [pp1x(:);flipud(pp2x(:))];
ppy = [pp1y(:);flipud(pp2y(:))];
%*************************************

%find particles that pass though polygon and are at the surface
I = find(inpolygon(xp,yp,ppx,ppy)&(bf < .1));

%all particles in the polygon
II = find(inpolygon(xp,yp,ppx,ppy));

%load SPM structure
SPM = SPMload;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
bed = SPM.data.bed;

Ib = find(inpolygon(Xc(:),Yc(:),ppx,ppy));
ys = 4000:50:7000;
for i=2:length(ys),
    Is = find((Yc(Ib) > ys(i-1))&(Yc(Ib)<ys(i)));
    bs(i-1) = min(bed(Ib(Is)));
end;
ys = .5*(ys(1:(end-1))+ys(2:end));


%[ys,Is] = sort(Yc(Ib));
%beds = bed(Ib);

set(gca,'position',[.1,.1,.8,.8]);
hold on; box on;
colormap('jet');

%grid structure
miny = 4000;
maxy = 7000;
minz = 4400;
maxz = 4800;
Nz = 100;
Ny = 100;
zint = linspace(minz,maxz,Nz); dz = zint(2)-zint(1);
yint = linspace(miny,maxy,Ny); dy = yint(2)-yint(1);
[ybin,zbin]=meshgrid(yint,zint);
Pnumber = zeros(size(ybin));
Psedi = zeros(size(ybin));
PN10 = zeros(size(ybin));
Page = zeros(size(ybin));
Pdist = zeros(size(ybin));
dzi = zint(2)-zint(1);
for i=1:length(II),

    jj = ceil((yp(II(i))-miny)/dy);
    ii = ceil((zp(II(i))-minz)/dz);
    
    if ((jj > 0)&(jj <= Ny)&(ii > 0)&(ii < Nz)), 
        Pnumber(ii,jj) = Pnumber(ii,jj) + 1;
        Psedi(ii,jj) = Psedi(ii,jj) + sedi(II(i));
        PN10(ii,jj) = PN10(ii,jj) + N10(II(i));
        Page(ii,jj) = Page(ii,jj) + age(II(i));
        Pdist(ii,jj) = Pdist(ii,jj) + dl(II(i));
    end;
    
end;


PN10 = PN10./(Pnumber+1e-6);
Page = Page./(Pnumber+1e-6);
Pcon = Psedi/dz;
Pdist = Pdist./(Pnumber+1e-6);

map = colormap;
map(1,:) = [1,1,1];
colormap(map);

subplot(2,2,1); hold on;
pcolor(yint,zint,PN10);
%contourf(ybin,zbin,PN10,20,'linestyle','none');
line(ys,bs,'color','k','linewidth',1);
shading flat; 
%caxis([20,200]*1e3);
cb = colorbar;
axis([4000,7000,4400,4800]);
xlabel('distance (m)');
ylabel('elevation (m)');
title('N10');

subplot(2,2,2);
pcolor(yint,zint,Page);
%contourf(ybin,zbin,PN10,20,'linestyle','none');
line(ys,bs,'color','k','linewidth',1);
shading flat; 
%caxis([20,50]*1e3);
cb = colorbar;
axis([4000,7000,4400,4800]);
xlabel('distance (m)');
ylabel('elevation (m)');
title('Age');

subplot(2,2,3);
pcolor(yint,zint,Pdist);
%contourf(ybin,zbin,PN10,20,'linestyle','none');
shading flat; 
line(ys,bs,'color','k','linewidth',1);
%caxis([20,50]*1e3);
cb = colorbar;
axis([4000,7000,4400,4800]);
xlabel('distance (m)');
ylabel('elevation (m)');
title('Travel distance');

subplot(2,2,4);
pcolor(yint,zint,Pcon);
%contourf(ybin,zbin,PN10,20,'linestyle','none');
shading flat; 
line(ys,bs,'color','k','linewidth',1);
caxis([0,1]);
cb = colorbar;
axis([4000,7000,4400,4800]);
xlabel('distance (m)');
ylabel('elevation (m)');
title('Debris concentration');

saveas(gcf,['./figs/psection',sprintf('%04d',fnr),'.fig'],'fig')


%Add sample positions
%load CNpoints.mat


if jpg,
    pause(.2);
    eval(['print -djpeg90 -r100 ./flic/flic',sprintf('%04d',fnr),'.jpg']);
end;


