function SPM = setup(L,H,nx,ny)

%
% default setup file for iSOSIA 3
%
% DLE 19/1 2015
%

%****** mesh data ******
SPM.mesh.L = L;          %mesh width
SPM.mesh.H = H;          %mesh height
SPM.mesh.nx = nx;        %number of cells
SPM.mesh.ny = ny;        %number of cells
SPM.mesh.ct = 0.1;       %time step scaling factor
SPM.mesh.maxtime = 1e3;  %length of simulation
SPM.mesh.filetime = 10;  %time between outputs
SPM.mesh.maxdt = 1.0;    %maximum time step
SPM.mesh.nmoni = 100;
SPM.mesh.gravity = 10;   %acceleration of gravity
SPM.mesh.periodic = 0;   %periodic boundary conditions
SPM.mesh.slidingmode = 2; %0:Weertman, 1:Budd, 2:Schoof
SPM.mesh.coldbased = 0;  %include coldbased ice 0/1
SPM.mesh.hydromode = 2;  %1:steady state, 2: transient
SPM.mesh.dofluvial = 0;
SPM.mesh.doice = 0;      %glacial flag
SPM.mesh.dosliding = 0;  %sliding flag
SPM.mesh.doperiglacial = 0; %periglacial flag
SPM.mesh.doglacialhydrology = 0; %glacial hydrology flag
SPM.mesh.doglacialerosion = 0; %glacial erosion flag
SPM.mesh.dohillslope = 0; %hillslope processes
SPM.mesh.dohillslopeerosion = 0; %sediment production on hillslopes
SPM.mesh.dolandslide = 0;  %landslide flag
SPM.mesh.doisostasy = 0; %isostasy flag
SPM.mesh.doavalance = 0; %avalance flag
SPM.mesh.docelldata = 0; %write time series from selected cells
SPM.mesh.dosediment = 0; %include layer of sediment
SPM.mesh.doglacialsedi = 0; %advect sediment volumes in glaciers
SPM.mesh.dodebrisablation = 0; %reduce surface ablation rates due
                               %to debris cover
SPM.mesh.maxb = 0.75;     %maximum bed gradient
SPM.mesh.maxs = 0.75;     %maximum ice surface gradient


%compute some extra parameters
SPM.mesh.dx = SPM.mesh.L/SPM.mesh.nx;
SPM.mesh.dy = SPM.mesh.H/SPM.mesh.ny;
SPM.mesh.hmin = min(SPM.mesh.dx,SPM.mesh.dy);


%****** iprop data ******
SPM.iprop.gamma = 1.0e-16*(910.0*9.81)^3.0; %Flow constant
SPM.iprop.gamma0 = 1.0e-12;  %Viscosity regelation
SPM.iprop.Cs = 1;
SPM.iprop.latentheat = 334e3;%latent heat for ice
SPM.iprop.ki = 2.14;         %Heat conductivity of ice
SPM.iprop.rho = 980;         %Density of ice
SPM.iprop.cp = 2000;         %Heat capacity of ice
SPM.iprop.ifac = 0.01;       %Relaxation constant (velocity)
SPM.iprop.sfac = 0.05;       %Relaxation constant (stress)
SPM.iprop.vbfac = 0.01;      %Relaxation constant (sliding)
SPM.iprop.maxitt_v = 500;    %Max iteration (velocity)
SPM.iprop.maxitt_s = 100;    %Max iteration (stress)
SPM.iprop.C = 0.42;          %Schoof sliding parameter
SPM.iprop.L0 = 4.0e-4;       %Schoof sliding parameter
SPM.iprop.minbeta = 0.02;    %Schoof sliding parameter
SPM.iprop.maxsliding = 300;  %Schoof sliding parameter
SPM.iprop.maxdeformation = 500; %max def. velocity
SPM.iprop.mf = 2.0;          %bedrock fracture density for
                             %quarrying
SPM.iprop.qfac = 1.0;        %pre-multiplication factor for quarrying
SPM.iprop.Ka = 1e-6;         %bedrock abrasion rate
SPM.iprop.ap = 2.0;          %abrasion sliding power

%***** mprop data *******
SPM.mprop.mtype = 2;
SPM.mprop.avaslope = 0.25; %max stable slope for avalances
SPM.mprop.avacurv = -0.1;  %min stable curvature for avalances (wind)
SPM.mprop.maxslope = 0.25; %max slope of avalanced snow
SPM.mprop.lrate = 0.006;   %lapse rate
SPM.mprop.Tsl = 0.0;       %temperature of slow line
SPM.mprop.maxacc = 2.0;    %max accumulation
SPM.mprop.maxabla = 10.0;  %max ablation rate
SPM.mprop.accgrad = 1.0;   %accumulation slope
SPM.mprop.ablgrad = 0.5;   %ablation slope
SPM.mprop.qb = 0.045;      %crustal heat flux
SPM.mprop.sedifrac = 1e-3; %debris fraction of accumulated ice
SPM.mprop.Ldebris = 0.72;  %ablation reduction beneath debris cover
                           %factor = exp(-S/Ldebris)
SPM.mprop.Temp = [0;10];   %t;T data
SPM.mprop.Mrate = [5000,6000;-1,1];   %h;Mrate data

%****** hwprop data ******
SPM.hwprop.Kgw = 1.0e3;    %water conductivity
SPM.hwprop.a2w = 0.01;     %fraction of ablation to glacial water 
SPM.hwprop.po = 0.01;      %ice porosity
SPM.hwprop.tscale = 1;     %time scaling < 1 slows system and
                           %increase timesteps

%**** hillslope properties **
SPM.hprop.Ks = 0.01;       %hillslope diffusivity
SPM.hprop.sc = 0.7;        %critical slope
SPM.hprop.Ke = 0.01;       %Hillslope erosion coefficient
SPM.hprop.Ls = 1.0;        %sediment thickness cover factor
SPM.hprop.gamma = 1;       %hillslope landslide strength
SPM.hprop.Nc = 20;         %number of random landslide cells

%**** fluvial properties ***
SPM.fprop.pr = 1.0;        %precipitation rate m/yr
SPM.fprop.rho_s = 2900;    %density of sediment grains
SPM.fprop.Dg = 0.01;       %average grain size
SPM.fprop.kw = 10.0;       %channel width scaling factor
SPM.fprop.tau_c = 0.01;    %critical shear stress
SPM.fprop.Kt = 11.1;       %transport capacity scaling
SPM.fprop.Ke = 1e-3;       %bedrock erosion rate


%****** grid data ******
SPM.data.bed = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.ice = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.sedi = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.maxacc = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.mrate = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.srate = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.fixflag = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.phi = 30*ones(SPM.mesh.ny,SPM.mesh.nx);

for i=1:20, SPM.data.Vs{i} = zeros(SPM.mesh.ny,SPM.mesh.nx); end;

%******** periglacial data **********
%if SPM.mesh.doperiglacial,
  load ./input/periglacial_input.mat;
  SPM.pgprop.nHs = length(Hsv);
  SPM.pgprop.nT0 = length(T0v);
  SPM.pgprop.Hsv = Hsv;
  SPM.pgprop.T0v = T0v;
  SPM.pgprop.Ci = Ci;
  SPM.pgprop.Tr = Tr;
  SPM.pgprop.rho_b = 2900.0;
  SPM.pgprop.rho_s = 2000.0;
  SPM.pgprop.Ke = 1.0e-3;
  SPM.pgprop.Kt = 1.0;         %Transport scaling
  SPM.pgprop.maxsedi = 10.0;   %Max sediment for periglacial erosion
  SPM.pgprop.maxice = 10.0;    %Max ice for periglacial erosion
  SPM.pgprop.minslope = 0.5;   %minimum slope for sediment free cells
  SPM.pgprop.minsedi = 1.0;    %expected minimum sediment cover on
                               %slopes smaller than minslope
  %end;



