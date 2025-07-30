function showhydro()

close all;
set(gcf,'units','centimeters','position',[10,10,25,40]);

%load data
tdata = load('./output/tseries.dat');
  
subplot(6,1,1);
hold on; box on;
plot(tdata(:,1),tdata(:,29),'-k');
ylabel('Melt input rate');
hold off;

subplot(6,1,2);
hold on; box on;
plot(tdata(:,1),tdata(:,26),'-b');
ylabel('hw (m)');
hold off;

subplot(6,1,3);
hold on; box on;
plot(tdata(:,1),tdata(:,27),'-r');
ylabel('dhw/dt (m)');
hold off;

subplot(6,1,4);
hold on; box on;
plot(tdata(:,1),tdata(:,30),'-k');
ylabel('sliding (m/yr)');
hold off;

subplot(6,1,5);
hold on; box on;
plot(tdata(:,1),tdata(:,25),'-k');
ylabel('N');
hold off;

subplot(6,1,6);
hold on; box on;
plot(tdata(:,1),tdata(:,9),'-k');
ylabel('Quarrying');
hold off;

save tdata.mat tdata 