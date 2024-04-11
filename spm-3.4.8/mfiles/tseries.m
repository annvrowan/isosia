function tseries()

close all
set(gcf,'units','centimeters','paperunits','centimeters');
set(gcf,'position',[5,15,20,20],'paperposition',[2,15,20,20]);

%load SPMstructure
SPM = SPMload;

%load data
tdata = load('./output/tseries.dat');

%sealevel temperature
subplot(2,1,1); hold on; box on;
line(tdata(:,1),tdata(:,10)-tdata(1,10),'color','k','linewidth',2);
axis([0,3000,-1,1]);
xlabel('time (yr)');
ylabel('dT (degrees)');

%Total ice volume
subplot(2,1,2); hold on; box on;
line(tdata(:,1),tdata(:,3)*SPM.mesh.dx*SPM.mesh.dy*SPM.mesh.nx*SPM.mesh.ny*1e-9,'color','k','linewidth',2);
axis([0,3000,.8,2]);
xlabel('time (yr)');
ylabel('Ice volume (km3)');

saveas(gcf,['./tseries.fig'],'fig')