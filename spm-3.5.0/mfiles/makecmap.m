function makecmap()

xc = linspace(0,1,256);

%topography part
%xt = linspace(1,0,6);
%rt = [225,200,195,154,148,245]/255;
%gt = [230,210,200,163,144,245]/255;
%bt = [215,190,180,144,133,245]/255;

%xt = linspace(0,1,4);
xt = [0.0,0.2,0.4,0.7,0.8,1.0];
rt = [0.6,0.4,0.3,0.8,0.6,0.6];
gt = [0.8,0.5,0.4,0.7,0.5,0.3];
bt = [0.5,0.3,0.2,0.6,0.5,0.2];

mapt = [interp1(xt,rt,xc)',...
        interp1(xt,gt,xc)',...
        interp1(xt,bt,xc)'];

%ice part
xi = linspace(0,1,2);
ri = [.9,.8];
gi = [.9,.8];
bi = [.99,.9];

mapi = [interp1(xi,ri,xc)',...
        interp1(xi,gi,xc)',...
        interp1(xi,bi,xc)'];

map = [mapt;mapi];


%sediment part
xi = linspace(0,1,2);
ri = [.9,.5];
gi = [.6,.3];
bi = [.2,.1];

maps = [interp1(xi,ri,xc)',...
        interp1(xi,gi,xc)',...
        interp1(xi,bi,xc)'];

%water part
xi = linspace(0,1,2);
ri = [0,0];
gi = [0,0];
bi = [.8,.5];

mapw = [interp1(xi,ri,xc)',...
        interp1(xi,gi,xc)',...
        interp1(xi,bi,xc)'];


%data part
colormap('default');
dmap = colormap;

ri = dmap(:,1);
gi = dmap(:,2);
bi = dmap(:,3);
xi = linspace(0,1,length(ri));

mapd = [interp1(xi,ri,xc)',...
        interp1(xi,gi,xc)',...
        interp1(xi,bi,xc)'];

map = [mapt;mapi;maps;mapw;mapd];

close all;
colormap(map);
colorbar

save cmap.mat map;
