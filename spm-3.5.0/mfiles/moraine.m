function zp=moraine(fnr,jpg)

close all

pcol1 = [.8,.4,.4];
pcol2 = [.5,.5,.5];

show('file',fnr,'withice','view',[0,90],'nocbar');
%show('file',fnr,'view',[0,90],'nocbar');
hold on;
set(gcf,'units','normalized','position',[.1,.2,.6,.6]);
%set(gcf,'paperunits','centimeters','paperposition',[2,5,20,16]);
set(gca,'position',[.05,.1,.4,.8]);
set(gca,'xlim',[800,7000]);
set(gca,'ylim',[900,10500]);

%load moraine positions
load mmoraine.mat;

%load particles data
loadparticles;


SPM = SPMload;

if 1,
    x0 = lx(1);
    y0 = ly(1);
    dw = 100; %width of 
    
    %**** form moraine polygon **********
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

    %line(ppx,ppy,8e3*ones(size(ppx)),'color','r');

    %find particles that pass though polygon and are at the surface
    I = find(inpolygon(xp,yp,ppx,ppy)&(bf < .1));
    II = find(inpolygon(xp,yp,ppx,ppy));
    
end;

load sourcearea.mat

nj(1) = round(sx(1)/SPM.mesh.dx);
nj(2) = round(sx(2)/SPM.mesh.dx);
ni(1) = round(sy(4)/SPM.mesh.dy);    
ni(2) = round(sy(1)/SPM.mesh.dy);    


plot3(xp(:),yp(:),6000*ones(size(xp(:))),'.','color',pcol2,'markersize',2);
%I = find((bpi>ni(1))&(bpi<ni(2))&(bpj>nj(1))&(bpj<nj(2))&(bf<0.1)); 
plot3(xp(I),yp(I),6000*ones(size(xp(I))),'.','color',pcol1,'markersize',2);
%plot3(xp(I),yp(I),zp(I),'.','color',pcol1,'markersize',.1);

%Add sample positions
load CNpoints.mat
hold on; 
zs = 7000*ones(size(xs));
for i=1:length(xs),
    plot3(xs,ys,zs,'db','markersize',10);
end;

%plot morraine 
axes('position',[.5,.5,.4,.4]);
hold on; box on;
plot(yp(I),1e-3*N10(I),'.','color',pcol1,'markersize',5);
axis([4000,7000,20,70]);
xlabel('distance (m)');
ylabel('N10 (x1000)');

errorbar(ys,1e-3*N10s,1e-3*dN10s,'.b');
plot(ys,1e-3*N10s,'db','markersize',10);



axes('position',[.5,.05,.4,.4]);
hold on; box on;
%plot(yp(II),zp(II),'.','color',pcol2,'markersize',2);

N10min = 20e3;
N10max = 50e3;
C = (N10(II)-N10min)/(N10max-N10min)+4; caxis([0,5]);
scatter(yp(II),zp(II),1*ones(size(yp(II))),C(:));
%plot(yp(I),zp(I),'.','color',pcol1,'markersize',2);
%cb = colorbar;
%set(cb,'ylim',[4,5]);
%axis([4000,7000,20,70]);
xlabel('distance (m)');
ylabel('elevation (m)');


if jpg,
    pause(.2);
    
    eval(['print -djpeg90 -r100 ./flic/flic',sprintf('%04d',fnr),'.jpg']);
    
end;


saveas(gcf,['./figs/moraine',sprintf('%04d',fnr),'.fig'],'fig')



