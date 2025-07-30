function ehist3(fnr)

close all;

set(gcf,'units','centimeters','paperunits','centimeters');
set(gcf,'position',[10,10,24,24],'paperposition',[8,15,8,8]);
set(gcf,'color','w');

SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;

ibed = SPM.data.bed;

load burial.mat;

Nbin = 35;
ebin = linspace(-1500,2300,Nbin+1); 
hbin = ebin;%linspace(-3000,1500,Nhypso+1); 
timefac = 700e3/84000;
filetime = 100;

Nf = length(ex(1,:));
erosion = zeros(Nbin,Nf);
erates = zeros(Nbin,Nf);
elevation = zeros(Nbin,Nf);
hypsoi = zeros(Nbin,1);
hypso = zeros(Nbin,1);
time = [0:(Nf-1)]*filetime;

iab = zeros(size(ibed));
fni = 0;

%colormap;
colormap('default');
map = colormap;
xx = linspace(0,1,length(map(:,1)));
xp = linspace(0,1,Nbin);
rp = interp1(xx,map(:,1)',xp);
gp = interp1(xx,map(:,2)',xp);
bp = interp1(xx,map(:,3)',xp);

loaddata; fbed = bed;
minbed = min(bed(:))
maxbed = max(bed(:))



%hypso - initial
for j = 1:Nbin,
    I = find((ibed(:) > hbin(j))&(ibed(:) <= hbin(j+1))&(xc(:) > 52e3)); 
    hypsoi(j) = length(I);
end;
hypsoi = hypsoi/(nx*ny);

%hypso - end

for j = 1:Nbin,
    I = find((fbed(:) > hbin(j))&(fbed(:) <= hbin(j+1))&(xc(:) > 52e3)); 
    hypso(j) = length(I);
end;
hypso = hypso/(nx*ny);

%erosion and erosion rates
for jj = 1:Nbin,
    I = find((fbed(:) > ebin(jj))&(fbed(:) <= ebin(jj+1))&(xc(:) > 52e3)); 
    if ~isempty(I),
        erosion(jj,:) = mean(ex(I,:),1); 
        erates(jj,2:end) = mean(diff(ex(I,:),1,2),1); 
        elevation(jj,:) = mean(topohist(I,:),1);
    else,
        erosion(jj,:) = 0; 
        erates(jj,2:end) = 0; 
        elevation(jj,:) = 0;
    end;
end;
erates = erates/timefac;


figure
hold on; box on;
for j= 2:Nbin,
    patch([0,-hypsoi(j),-hypsoi(j),0],[hbin(j-1),hbin(j-1),hbin(j),hbin(j)],[rp(j),gp(j),bp(j)]);
    patch([0,hypso(j),hypso(j),0],[hbin(j-1),hbin(j-1),hbin(j),hbin(j)],[rp(j),gp(j),bp(j)]);
end;





figure
hold on; box on;
for j = 1:Nbin,
    line(time,erosion(j,:),'color',[rp(j),gp(j),bp(j)]);
end;

return


%erosion rates
axes('position',[.05,0.05,0.9,0.9]);
hold on; box on;
set(gca,'color','w');
%line(erates(2:(Nbin-1),1),elevation(2:(Nbin-1),end),'color','k')
%line(erates(2:(Nbin-1),end),elevation(2:(Nbin-1),end),'color','k')
I = find(erates(:) < 1e-6); erates(I) = 1e-6;
for j = 2:1:(Nbin-1),
    line(erates(j,:),.5*(elevation(j,:)+elevation(j-1,:)),'color',[rp(j),gp(j),bp(j)],'linewidth',.5);
    plot(erates(j,1),.5*(elevation(j,1)+elevation(j-1,1)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[rp(j),gp(j),bp(j)],'linewidth',.5); 
    plot(erates(j,end),.5*(elevation(j,end)+elevation(j-1,end)),'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[rp(j),gp(j),bp(j)],'linewidth',.5); 
end;
%axis([1e-6,1e-2,-1500,1500]);
set(gca,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2]);
set(gca,'ytick',[-1500:500:1000]);
%set(gca,'xticklabel',[],'yticklabel',[]);
set(gca,'xscale','log');


return

save ehist.mat erates elevation erosion ebin

print -depsc -loose flic/erates.eps
print -djpeg90 -r300 flic/erates.jpg


