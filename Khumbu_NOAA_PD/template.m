function template()

% close all;

%Define mesh dimensions
L = 27.2e3; H = 13.8e3;        %Dimensions of domain in km
nx = 272; ny = 138;            %No. cells in each axis

%load default SPM structure
SPM = setup(L,H,nx,ny);
SPM.mesh.dosliding = 1;
SPM.mesh.slidingmode = 2; %0:Weertman, 1:Budd, 2:Schoof
SPM.mesh.doperiglacial = 0;
SPM.mesh.doglacialerosion = 0;
SPM.mesh.doavalance = 1;
SPM.mesh.doglacialhydrology = 0;
SPM.mesh.doglacialsedi = 1; %Sediment transport within the ice
SPM.mprop.sedifrac = 1e-2;  %debris fraction of accumulated ice
SPM.mesh.dosediment = 1;    %Sediment deposition at ice margins
SPM.mesh.dohillslope = 1;   %hillslope processes
SPM.mesh.nmoni = 10;
SPM.mesh.maxtime = 120;    %length of simulation
SPM.mesh.filetime = 10;     %time between outputs
SPM.mesh.maxdt = 1;         %maximum time step
SPM.mesh.ct = 0.1;
SPM.mesh.dodebrisablation = 1; %Adjust ablation rate in line with debris cover

%Hydrological switches
SPM.mesh.doice = 1;
SPM.data.fixflag(1,:) = 1;
SPM.data.fixflag(end,:) = 1;
SPM.data.fixflag(:,1) = 1;
SPM.data.fixflag(:,end) = 1;

%% Grid DEM for model domain
%Mesh the model domain
x = linspace(0,SPM.mesh.L,SPM.mesh.nx);
y = linspace(0,SPM.mesh.H,SPM.mesh.ny);
[X,Y] = meshgrid(x,y);
load 'khumbu_100m_dem.mat';
bed = DEM.Z;
load 'icediff_6.mat';
bed = bed+icediff_6;
SPM.data.bed = bed;

%% Define the extent of the catchment
py=5;
load 'khumbu_100m_domain.mat';
I = inpolygon(X,Y,px,py);
SPM.data.maxacc(I) = 10.0;


%% Mass Balance switches
SPM.mprop.mtype = 3;       %Type of MB calculation to use: 1 = vectors for MB-Z; 2 = Gradients for acc/abl; 3 = gridded EMB data
      
SPM.mprop.avaslope = 0.50; %max stable slope for avalances
SPM.mprop.maxslope = 0.25; %max slope of avalanced snow

SPM.mprop.Ldebris = 0.8;

%Mtype = 3. Read in EMB data
load 'khumbu_100_20952100_mb_noaa_rcp45.mat';
meltrate(isnan(meltrate))=0;
SPM.data.mrate = meltrate;


%% Load in existing ice thickness
%load data
    fid = fopen('PD_icemass_o20.dat');
    xc = fread(fid,[ny,nx],'double');
    yc = fread(fid,[ny,nx],'double');
    bed = fread(fid,[ny,nx],'double');
    ice = fread(fid,[ny,nx],'double');
    bslope = fread(fid,[ny,nx],'double');
    tslope = fread(fid,[ny,nx],'double');
    te2 = fread(fid,[ny,nx],'double');
    frost = fread(fid,[ny,nx],'double');
    sediment = fread(fid,[ny,nx],'double');
    tn = fread(fid,[ny,nx],'double');
    te = fread(fid,[ny,nx],'double');
    ts = fread(fid,[ny,nx],'double');
    sliding = fread(fid,[ny,nx],'double');
    deformation = fread(fid,[ny,nx],'double');
    pw = fread(fid,[ny,nx],'double');
    bmelt = fread(fid,[ny,nx],'double');
    quarrying = fread(fid,[ny,nx],'double');
    accrate = fread(fid,[ny,nx],'double');
    meltrate = fread(fid,[ny,nx],'double');
    frostrate = fread(fid,[ny,nx],'double');
    fluvial = fread(fid,[ny,nx],'double');
    landslide = fread(fid,[ny,nx],'double');
    abrasion = fread(fid,[ny,nx],'double');
    isostasy = fread(fid,[ny,nx],'double');
    lee = fread(fid,[ny,nx],'double');
    Ts = fread(fid,[ny,nx],'double');
    hillslope = fread(fid,[ny,nx],'double');
    Ta = fread(fid,[ny,nx],'double');
    Tb = fread(fid,[ny,nx],'double');
    sfac = fread(fid,[ny,nx],'double');
    margin = fread(fid,[ny,nx],'int');
    for i=1:20,
        Vs{i} = fread(fid,[ny,nx],'double');
    end
    fclose(fid);

SPM.data.ice = ice;
SPM.data.Vs = Vs;

%% Iceflow params
SPM.iprop.ifac = 0.002;      %Relaxation constant (velocity)
SPM.iprop.sfac = 0.02;       %Relaxation constant (stress)
SPM.iprop.maxitt_v = 500;    %Max iteration (velocity)
SPM.iprop.maxitt_s = 200;    %Max iteration (stress)
SPM.iprop.gamma = 1.0e-16*(910.0*9.81)^3.0; %Flow constant
SPM.iprop.C = 0.42;          %Schoof sliding parameter
SPM.iprop.L0 = 4.0e-4;       %Schoof sliding parameter

%**** hillslope properties **
SPM.hprop.Ks = 0.1;       %hillslope diffusivity
SPM.hprop.sc = 0.7        %critical slope

surf(X,Y,bed); colorbar, axis tight

write(SPM);


