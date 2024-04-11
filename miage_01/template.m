function template()

close all;

%mesh dimensions
L = 14.9e3; H = 10.9e3;        %Dimensions of domain in km
nx = 300; ny = 220;            %No. cells in each axis
SPM = setup(L,H,nx,ny);

%****** mesh data ******
SPM.mesh.doice = 1;
SPM.mesh.dosliding = 1;  %sliding flag
SPM.mesh.doavalance = 1; %avalance flag
SPM.mesh.slidingmode = 2; %0:Weertman, 1:Budd, 2:Schoof
SPM.mprop.dhice = 0;     %ice contribution to elevation used to compute mass balance (0 = bed, 1 = ice surface)

%Switches for ice and sediment functions
SPM.mesh.dohillslope = 1; %hillslope processes
SPM.mesh.dohillslopeerosion = 1; %bedrock erosion on hillslopes
SPM.mesh.doweathering = 0; %sediment production on hillslopes
SPM.mesh.dosediment = 1; %include layer of sediment
SPM.mesh.doglacialerosion = 1; %glacial erosion flag
SPM.iprop.Ka = 1e-6;         %bedrock abrasion rate
SPM.mesh.doglacialsedi = 1; %advect sediment volumes in glaciers
SPM.mesh.doparticles = 1; %use lagrangian particles to advect sediment
SPM.mesh.dodebrisablation = 1; %reduce surface ablation rates due to debris cover
SPM.mesh.dodeposit = 1; %allows particles to form sediment when killed

%****** iprop data ******
SPM.iprop.ifac = 0.002;       %Relaxation constant (velocity)
SPM.iprop.sfac = 0.02;       %Relaxation constant (stress)
SPM.iprop.vbfac = 0.001;      %Relaxation constant (sliding)
SPM.iprop.maxitt_v = 500;    %Max iteration (velocity)
SPM.iprop.maxitt_s = 500;    %Max iteration (stress)
SPM.iprop.C = 0.42;          %Schoof sliding parameter
SPM.iprop.L0 = 4.0e-4;       %Schoof sliding parameter
SPM.iprop.minbeta = 0.005;    %Schoof sliding parameter
SPM.iprop.maxsliding = 200;  %Schoof sliding parameter
SPM.iprop.maxdeformation = 200; %max def. velocity
SPM.iprop.Ka = 1e-6;         %bedrock abrasion rate
SPM.iprop.ap = 2.0;          %abrasion sliding power

%***** Sediment particles ************** 
SPM.parprop.npmax = 1e9; %max number of particles
SPM.parprop.minice = 10.0; %minimum ice thickness for particle formation
SPM.parprop.minsedi = 0.05; %minimum sediment thickness for particle formation
SPM.parprop.maxsedi = 0.2; %max sediment per particle
SPM.parprop.maxpm = 100; %max particle number in margin cell
SPM.parprop.maxp = 200; %max particle number in ordinary cell
SPM.parprop.minpm = 0; %minimum particle number at margin 
SPM.parprop.minsedim = 0; %mimum sediment pickup at margin
                          %The two latter parameters are used to
                          %form sediment along margins without
                          %production of sediment in the landscape


%% Set sim time, define mass balance and read in domain
SPM.mesh.maxtime = 14.1e3;  %length of simulation
SPM.mesh.filetime = 200;  %time between outputs

%******* mprop *********
SPM.mprop.mtype = 2; %1 = elevation dependent, 2 = temperature dependent, 3 = EMB?
SPM.mprop.avaslope = 0.5; %max stable slope for avalances
SPM.mprop.avacurv =  -0.1;  %min stable curvature for avalances (wind)
SPM.mprop.lrate = 0.006;   %lapse rate
SPM.mprop.maxacc = 2.0;  %maximum accumulation rate
SPM.mprop.mPDD = 2e-3;     %Positive degree day factor (see massbalance.m)
SPM.mprop.dTa = 6;        %Annual temperature vatiation (see massbalance.m)
SPM.mprop.Ldebris = 0.2;   %ablation reduction beneath debris cover factor = exp(-S/Ldebris)
precipfactor = 1;

%**** hillslope properties **
SPM.hprop.Ks = 0.01;       %hillslope diffusivity
SPM.hprop.sc = 0.5;        %critical slope
SPM.hprop.Ke = 0.001;       %Hillslope erosion coefficient
SPM.hprop.Kw = 1e-2;       %Hillslope weathering coefficient w=Kw*exp(-Hs/Ls)
SPM.hprop.Ls = 1.0;        %sediment thickness cover factor
SPM.hprop.gamma = 1;       %hillslope landslide strength
SPM.hprop.Nc = 20;         %number of random landslide cells 
SPM.hprop.maxsedi = 10;    %max sediment thickness on hillslopes

%Temperature forcing from temp12k 30-60N (Kaufmann et al., 2020)
load 'temp12k_spm.mat';
t2 = temp12k_spm.t + 1000;
T2 = temp12k_spm.T_3060N_Ann + 12;
t1 = [0 1000];
T1 = [9.4 9.4];
t3 = [t2(1,end),t2(1,end)+1000];
T3 = [T2(1,end)+1,T2(1,end)+1];
t = cat(2,t1,t2,t3);
T = cat(2,T1,T2,T3);

%Set into SPM
SPM.mprop.Temp = [t;T];    %t;T data


%********** Load bed topo ******************
x = linspace(0,SPM.mesh.L,SPM.mesh.nx);
y = linspace(0,SPM.mesh.H,SPM.mesh.ny);
[X,Y] = meshgrid(x,y);

%Load DEM and resample bed
load miage_dem_huss.mat;
bed = interp2(DEM.X,DEM.Y,DEM.Z,X,Y,'linear',1000);
SPM.data.bed = bed;


SPM.data.fixflag(:,1) = 1;
%Define domain with include
py=5;
load 'miage_domain.mat';
SPM.data.include = inpolygon(X,Y,px,py);

%% Look at set up for this simulation
figure(1)
surf(X,Y,bed);
figure(2)
plot(t,T)

%mass balance function (positive degree day)
SPM = massbalance(SPM,1);

write(SPM);