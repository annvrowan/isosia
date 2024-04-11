function SPM=massbalance(SPM,fg)

%computes and plots massbalance models%
% SPM is SPM structure
% fg = 0 : no graphics
% fg = 1 : make graphics
%

%close all;

%load SPM data
mPDD = SPM.mprop.mPDD; %Positive degree day melt rate constant
dTa = SPM.mprop.dTa; %annual temperature variation amplitude
precip = 1.0;%acc rates are nomalised to precipitation rate = 1 m/y
             %local acc rates must therefore be scaled by
             %SPM.data.precip inside isosia
Tsl = SPM.mprop.Temp(2,1); %Sea level temperture
lrate = SPM.mprop.lrate; %lapse rate

%input 
nt = 365; %days per year
nd = 200; %number of temperature points

%Temperature range
T0 = linspace(-2*dTa,2*dTa,nd);

%Time vector - one year
time = [1:365];

%compute temperature matrix
Temp = T0(:)*ones(1,nt) + ones(nd,1)*dTa*sin(2*pi*time(:)'/365);

%Compute positive degree day index and number of days with frost
PDD = zeros(nd,1);
nFd = zeros(nd,1);
for i=1:nd,
    I = find(Temp(i,:) > 0);
    PDD(i) = sum(Temp(i,I));
    nFd(i) = nt - length(I);
end;

%melt rate for varying T0
melt = mPDD*PDD;

%accumulation rate - number of days with frost times average precip rate
macc = precip*nFd/nt;

%massbalance
mrate = macc - melt;

%Associated elevation range - using first temperature input
h = (Tsl - T0)/lrate;

%save massbalance data to SPM structure
SPM.mprop.Mrate_T = [T0(:)';macc(:)';melt(:)'];

%Print out location of ELA (AVR)
ela_mrate = round(mrate,1);
ela_i = find(ela_mrate == 0);
ela = mean(h(ela_i))

figure
subplot(1,2,1);
hold on; box on; grid on;
li(1)=line(T0,-melt,'color','r');
li(2)=line(T0,macc,'color','b');
li(3)=line(T0,mrate,'color','k');
legend(li,'Melt','Accumulation','Mass balance','Location','SouthWest');
xlabel('Temperature (C)');
ylabel('Mass balance (m/y)');

subplot(1,2,2);
hold on; box on; grid on;
li(1)=line(-melt,h,'color','r');
li(2)=line(macc,h,'color','b');
li(3)=line(mrate,h,'color','k');
legend(li,'Melt','Accumulation','Mass balance','Location','NorthWest');
ylabel('Elevation (m)');
xlabel('Melt rate (m/y)');

