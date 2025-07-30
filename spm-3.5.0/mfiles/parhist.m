function zp=parhist(fnr,jpg)


%prepare flowlines
prepare_flowlines(fnr);
flowlines;

close all
set(gcf,'units','normalized','position',[.1,.2,.6,.6]);

pcol1 = [.8,.4,.4];
pcol2 = [.5,.5,.5];

%load moraine positions
load mmoraine.mat;

%load particles data
loadparticles;

SPM = SPMload;

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

load flines.mat;

%loop flowlines
for i = 1:length(flines),
    I = find(inpolygon(flines{i}.x,flines{i}.y,ppx,ppy)),
    if any(flines{i}.bf(I) < 0.05),
        line(flines{i}.y,flines{i}.z,'color','k','linewidth',1);
    end;
end;


%line(ppx,ppy,8e3*ones(size(ppx)),'color','r');

%find particles that pass though polygon and are at the surface
I = find(inpolygon(xp,yp,ppx,ppy)&(bf < .1));
II = find(inpolygon(xp,yp,ppx,ppy));

%Add sample positions
load CNpoints.mat

set(gca,'position',[.1,.1,.8,.8]);
hold on; box on;
colormap('jet');

C = N10(II); caxis([20,50]*1e3);
scatter(yp(II),zp(II),1*ones(size(yp(II))),C(:));
%plot(yp(I),zp(I),'.','color',pcol1,'markersize',2);
cb = colorbar;
%set(cb,'ylim',[4,5]);
axis([3500,8000,4350,4900]);
xlabel('distance (m)');
ylabel('elevation (m)');


if jpg,
    pause(.2);
    eval(['print -djpeg90 -r100 ./flic/flic',sprintf('%04d',fnr),'.jpg']);
end;


