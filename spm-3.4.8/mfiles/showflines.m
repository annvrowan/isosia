function showflines()

close all

load flines.mat;
load mmoraine.mat;
x0 = lx(1);
y0 = ly(1);
dw = 50;

% form polygon
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

%find flowlines passing through polygon
mor = [];
nhit = 0;
for i=1:length(flines), %loop all flowlines
    x = flines{i}.x(:);
    y = flines{i}.y(:);
    z = flines{i}.z(:);
    bz = flines{i}.bz(:);
    time = flines{i}.time(:);
    I = find(inpolygon(x,y,ppx,ppy)); %find part of flowline in polygon
    if (length(I) > 1)&(any(bz(I) < 10)), %only continue if any
                                          %part of the flowline is
                                          %near surface within polygon
        nhit = nhit + 1;
        mor(nhit) = i; %registre this flowline
        II = find((inpolygon(x,y,ppx,ppy)==1)&(bz < 10)); 
        flowlines{i}.mi = II;
    end;    
end;

figure(1); hold on; box on;
for i=1:length(flines),
    lx = flines{i}.x;
    ly = flines{i}.y;
    lz = flines{i}.z;
    plot3(lx,ly,lz);    
end;
view([40,20]);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

savefig(gcf,'flowlines');



figure(2); hold on; box on;
for i=1:length(flines),
    lx = flines{i}.x;
    ly = flines{i}.y;
    lz = flines{i}.z;
    line(lx,ly,'color','k','linewidth',.5);    
end;

for i=1:length(mor),
    lx = flines{mor(i)}.x;
    ly = flines{mor(i)}.y;
    lz = flines{mor(i)}.z;
    line(lx,ly,'color','b','linewidth',1)    
end;
line(ppx,ppy,'color','r','linewidth',2);
load CNpoints.mat
plot(xs,ys,'.m','markersize',20);
xlabel('x (m)');
ylabel('y (m)');
savefig(gcf,'selected_lines');


figure(3); hold on; box on;
for i=1:length(mor),
    lx = flines{mor(i)}.x;
    ly = flines{mor(i)}.y;
    lz = flines{mor(i)}.z;
    ltime = flines{mor(i)}.time;
    ldist = flines{mor(i)}.dist;
    lNBe = flines{mor(i)}.NBe;
    [mind,nmin] = min(sqrt((lx-x0).^2+(ly-y0).^2));
    ldist = ldist - ldist(nmin);
    mi = flowlines{mor(i)}.mi;
    line(ldist(mi),lNBe(mi),'color','b')    
    %line(ltime(mi),lNBe(mi),'color','b')    
end;
xlabel('Distance in polygon (m)');
ylabel('NBe atoms/gram');
savefig(gcf,'Cosmo');




