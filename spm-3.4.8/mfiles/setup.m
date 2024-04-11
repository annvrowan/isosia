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
SPM.mesh.gmode = 1;      %glacier mode, 0:SIA, 1:iSOSIA
SPM.mesh.slidingmode = 3; %0:Weertman, 1:Budd, 2:Schoof 3: cavity model
SPM.mesh.coldbased = 0;  %include coldbased ice 0/1
SPM.mesh.hydromode = 3;  %1:steady state, 2: transient, 3: cavity/channel
SPM.mesh.dofluvial = 0;
SPM.mesh.doice = 0;      %glacial flag
SPM.mesh.dosliding = 0;  %sliding flag
SPM.mesh.doperiglacial = 0; %periglacial flag
SPM.mesh.doglacialhydrology = 0; %glacial hydrology flag
SPM.mesh.doglacialerosion = 0; %glacial erosion flag
SPM.mesh.dohillslope = 0; %hillslope processes
SPM.mesh.dohillslopeerosion = 0; %bedrock erosion on hillslopes
SPM.mesh.doweathering = 0; %sediment production on hillslopes
SPM.mesh.dolandslide = 0;  %landslide flag
SPM.mesh.doisostasy = 0; %isostasy flag
SPM.mesh.doavalance = 0; %avalance flag
SPM.mesh.doorographic = 0; %flag for orographic precipitation
SPM.mesh.docelldata = 0; %write time series from selected cells
SPM.mesh.dosediment = 0; %include layer of sediment
SPM.mesh.doglacialsedi = 0; %advect sediment volumes in glaciers
SPM.mesh.doparticles = 0; %use lagrangian particles to advect
                          %sediment
SPM.mesh.dodeposit = 0; %allows particles to form sediment when killed
SPM.mesh.dodebrisablation = 0; %reduce surface ablation rates due
                               %to debris cover
SPM.mesh.maxb = 10;     %maximum bed gradient
SPM.mesh.maxs = 10;     %maximum ice surface gradient

%compute some extra parameters
SPM.mesh.dx = SPM.mesh.L/SPM.mesh.nx;
SPM.mesh.dy = SPM.mesh.H/SPM.mesh.ny;
SPM.mesh.hmin = min(SPM.mesh.dx,SPM.mesh.dy);


%****** iprop data ******
SPM.iprop.gamma = 1.0e-16*(910.0*9.81)^3.0; %Flow constant
SPM.iprop.gamma0 = 1.0e-2;  %Viscosity regelation
SPM.iprop.Cs = 1e-3;
SPM.iprop.latentheat = 334e3;%latent heat for ice
SPM.iprop.ki = 2.14;         %Heat conductivity of ice
SPM.iprop.rho = 910;         %Density of ice
SPM.iprop.cp = 2000;         %Heat capacity of ice
SPM.iprop.ifac = 0.01;       %Relaxation constant (velocity)
SPM.iprop.sfac = 0.05;       %Relaxation constant (stress)
SPM.iprop.vbfac = 0.01;      %Relaxation constant (sliding)
SPM.iprop.maxitt_v = 500;    %Max iteration (velocity)
SPM.iprop.maxitt_s = 100;    %Max iteration (stress)
SPM.iprop.C = 0.42;          %Schoof sliding parameter
SPM.iprop.L0 = 4.0e-4;       %Schoof sliding parameter
SPM.iprop.minbeta = 0.02;    %Schoof sliding parameter
SPM.iprop.mf = 2.0;          %bedrock fracture density for
                             %quarrying
SPM.iprop.Kq = 1.0;        %pre-multiplication factor for quarrying
SPM.iprop.Ka = 1e-4;         %bedrock abrasion rate
SPM.iprop.ap = 2.0;          %abrasion sliding power
SPM.iprop.minefac = 1.0e-6;  %min abrasion constant
SPM.iprop.ksg = 1e-7;        %apparent conductivity of the array of
                             %debris to ice 

%***** mprop data *******
SPM.mprop.mtype = 2;
SPM.mprop.avaslope = 0.25; %max stable slope for avalances
SPM.mprop.avacurv = -0.1;  %min stable curvature for avalances
                           %(wind)
SPM.mprop.maxacc = 100.0;  %maximum accumulation rate
SPM.mprop.dhice = 1.0;     %ice contribution to elevation used to
                           %compute mass balance
SPM.mprop.lrate = 0.006;   %lapse rate
SPM.mprop.mPDD = 1e-3;     %Positive degree day factor (see massbalance.m)
SPM.mprop.dTa = 8.0;       %Annual temperature vatiation (see massbalance.m)
SPM.mprop.qb = 0.045;      %crustal heat flux
SPM.mprop.Ldebris = 0.72;  %ablation reduction beneath debris cover
                           %factor = exp(-S/Ldebris)
SPM.mprop.Temp = [0;10];   %t;T data - temperature at sealevel
SPM.mprop.Mrate_h = [5000,6000;-1,1];   %[h;Mrate] data - use with mprop=1
SPM.mprop.Mrate_T = [-10,10;1,1;1,1];   %[T;acc;melt] data - use
                                        %with mprop=2
SPM.mprop.dL = 10e3;       %smoothing length scale for orographic
                           %precip
SPM.mprop.precip = 1.0;    %precipitation rate
SPM.mprop.cp = 2.0e-3;     %precipitation constant
SPM.mprop.ch = 1.0e-6;     %orographic constant

%****** hwprop data ******
SPM.hwprop.a2w = 1;     %fraction of ablation to glacial water 
SPM.hwprop.tscale = 1;     %time scaling < 1 slows system and
                           %increase timesteps
SPM.hwprop.kh = 0.1506;    %run code htest.m to find these
SPM.hwprop.h0 = 0.1;  
SPM.hwprop.ds = 0.25;
SPM.hwprop.ss0 = 0.3;
SPM.hwprop.kmin = 1e-7;    %min transmicity for cavities 
SPM.hwprop.alpha = 0.7;    %cavity area scaling 
SPM.hwprop.Ls = 4;         %Cavity spacing
SPM.hwprop.Lc = 200;       %Channel spacing
SPM.hwprop.minqw = 1e-4;   %max water flux
SPM.hwprop.S0 = 0.1;       %min cavity length
SPM.hwprop.A0 = 0.1;       %min channel cross section
SPM.hwprop.dtw = 0.1;      %time step scaling
SPM.hwprop.minqw = 1.0e-4; %minimum water flux

%**** hillslope properties **
SPM.hprop.Ks = 0.01;       %hillslope diffusivity
SPM.hprop.sc = 0.7;        %critical slope
SPM.hprop.Ke = 0.01;       %Hillslope erosion coefficient
SPM.hprop.Kw = 1e-4;       %Hillslope weathering coefficient w=Kw*exp(-Hs/Ls)
SPM.hprop.Ls = 1.0;        %sediment thickness cover factor
SPM.hprop.gamma = 1;       %hillslope landslide strength
SPM.hprop.Nc = 20;         %number of random landslide cells 
SPM.hprop.maxsedi = 10;    %max sediment thickness on hillslopes

%**** fluvial properties ***
SPM.fprop.pr = 1.0;        %precipitation rate m/yr
SPM.fprop.rho_s = 2900;    %density of sediment grains
SPM.fprop.Dg = 0.01;       %average grain size
SPM.fprop.kw = 10.0;       %channel width scaling factor
SPM.fprop.tau_c = 0.01;    %critical shear stress
SPM.fprop.Kt = 11.1;       %transport capacity scaling
SPM.fprop.Ke = 1e-3;       %bedrock erosion rate

%***** Flexural isostasy ********
SPM.flexprop.Te = 10.0e3;  %Elastic thickness
SPM.flexprop.rho_r = 2800; %Density of eroded bedrock
SPM.flexprop.rho_s = 2400; %Density of sediment
SPM.flexprop.rho_a = 3250; %Density of asthenosphere
SPM.flexprop.ncall = 10;   %Call isostasy every ncall timestep
 

%***** Sediment particles ************** 
SPM.parprop.npmax = 1e5; %max number of particles
SPM.parprop.minice = 10.0; %minimum ice thickness for ice transport
SPM.parprop.minsedi = 0.05; %minimum sediment thickness for particle formation
SPM.parprop.maxsedi = 0.2; %max sediment per particle
SPM.parprop.maxpm = 100; %max particle number in margin cell
SPM.parprop.maxp = 200; %max particle number in ordinary cell
SPM.parprop.minpm = 0; %minimum particle number at margin 
SPM.parprop.minsedim = 0; %mimum sediment pickup at margin
                          %The two latter parameters are used to
                          %form sediment along margins without
                          %production of sediment in the landscape

%****** grid data ******
SPM.data.bed = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.ice = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.sedi = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.precip = ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.include = ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.mrate = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.srate = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.fixflag = zeros(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.phi = SPM.hprop.sc*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.CNprod = 4.0*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.atten = 150*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.Ka = SPM.iprop.Ka*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.Kq = SPM.iprop.Kq*ones(SPM.mesh.ny,SPM.mesh.nx);
SPM.data.Cs = SPM.iprop.Cs*ones(SPM.mesh.ny,SPM.mesh.nx);

%for i=1:20, SPM.data.Vs{i} = zeros(SPM.mesh.ny,SPM.mesh.nx); end;

%******** periglacial data **********
%if SPM.mesh.doperiglacial,
  load ./input/periglacial_input.mat;
  SPM.pgprop.nHs = length(Hsv);
  SPM.pgprop.nT0 = length(T0v);
  SPM.pgprop.Hsv = Hsv;
  SPM.pgprop.T0v = T0v;
  SPM.pgprop.Ci = Ci;
  SPM.pgprop.CiT = CiT;
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



