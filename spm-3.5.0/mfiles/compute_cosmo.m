function compute_cosmo(fnr)

close all;

SPM = SPMload;
L = SPM.mesh.L;
H = SPM.mesh.H;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

%load burial data
load burial.mat;
[ny,nx,Nt] = size(burial);

%************* Cosmo info **************

%Cosmogenic props
TBe = 1.378e6;
TAl = 0.705e6;
CNprop.lambda_Be = log(2)/TBe;
CNprop.lambda_Al = log(2)/TAl;
CNprop.rho = 2.65; %density
CNprop.PBe0 = 4.0; %nomalization production
CNprop.pratio = 7.1; %poduction ratio
CNprop.pr_fm_Be = 0.005; %fast muons (of total p)
CNprop.pr_fm_Al = 0.006;
CNprop.pr_nmc_Be = 0.015;
CNprop.pr_nmc_Al = 0.018;
CNprop.Lspal = 150; %g/cm^2
CNprop.Lnmc = 1500;
CNprop.Lfm = 4320;

%production parameters
P10tot = CNprop.PBe0;
P10spal = (1-CNprop.pr_fm_Be-CNprop.pr_nmc_Be)*P10tot;
P10fm = CNprop.pr_fm_Be*P10tot;
P10nmc = CNprop.pr_nmc_Be*P10tot;
P26tot = CNprop.pratio*CNprop.PBe0;
P26spal = (1-CNprop.pr_fm_Al-CNprop.pr_nmc_Al)*P26tot;
P26fm = CNprop.pr_fm_Al*P26tot;
P26nmc = CNprop.pr_nmc_Al*P26tot;

%steady state erosion rate
erate = 10e-6;            

ftime = 500;
time = [0:(fnr-1)]*ftime;
ntime = length(time);

TaBe = zeros(ny,nx);
NBe_M = zeros(ny,nx);
NAl_M = zeros(ny,nx);
noice = zeros(ny,nx);
h = waitbar(0,'Computing...');

for i=1:ny,
    for j=1:nx,
 
        
        %erosion history - increasing through time
        ehist = zeros(fnr,1);
        ehist(:) = ex(i,j,1:fnr);        
        
        %burial history  - decreasing in time and ending in 0
        epath = zeros(fnr,1);
        epath(:) = ehist(fnr) - ehist(:); 
        
        %ice cover history
        ihist = zeros(fnr,1);
        ihist(:) = icehist(i,j,1:fnr);
        if all(ihist(:) < 50), noice(i,j) = 1; end;
        
        %topography of surface through time
        thist = zeros(fnr,1);
        thist(:) = topohist(i,j,1:fnr);
        
        
        %production function
        pfac = ones(size(epath));
        I = find(ihist > 30); pfac(I) = 0; %ice cover
        I = find(thist < -30); pfac(I) = 0; %under sealevel
        
        %subplot(4,1,1); hold on; box on;
        %line(time,ehist,'color','k'); 
        %title('Ehist');
        
        %subplot(4,1,2); hold on; box on;
        %line(time,epath,'color','k'); 
        %title('Epath');
        
        %subplot(4,1,3); hold on; box on;
        %line(time,ihist,'color','k'); 
        %title('ihist');

        %subplot(4,1,4); hold on; box on;
        %line(time,thist,'color','k'); 
        %title('thist');

        %allocate
        N10 = zeros(ntime,1);
        N26 = zeros(ntime,1);
        
        %compute steady state starting concentration
        fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lspal;
        fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lspal;
        N10(1) = P10spal*exp(-CNprop.rho*100*epath(1)/CNprop.Lspal)/fBe;
        N26(1) = P26spal*exp(-CNprop.rho*100*epath(1)/CNprop.Lspal)/fAl;
        fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lnmc;
        fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lnmc;
        N10(1) = N10(1) + P10nmc*exp(-CNprop.rho*100*epath(1)/CNprop.Lnmc)/fBe;
        N26(1) = N26(1) + P26nmc*exp(-CNprop.rho*100*epath(1)/CNprop.Lnmc)/fAl;
        fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lfm;
        fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lfm;
        N10(1) = N10(1) + P10fm*exp(-CNprop.rho*100*epath(1)/CNprop.Lfm)/fBe;
        N26(1) = N26(1) + P26fm*exp(-CNprop.rho*100*epath(1)/CNprop.Lfm)/fAl;
        
        %integrate time
        for kk=2:ntime,
    
            dt = (time(kk)-time(kk-1));
            pf = .5*(pfac(kk)+pfac(kk-1));
            bz = .5*(epath(kk)+epath(kk-1));
            P10z = pf*(P10spal*exp(-CNprop.rho*100*bz/CNprop.Lspal)+P10nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P10fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
            P26z = pf*(P26spal*exp(-CNprop.rho*100*bz/CNprop.Lspal)+P26nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P26fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
            N10(kk) = N10(kk-1)*exp(-dt*CNprop.lambda_Be) + P10z*(1-exp(-dt*CNprop.lambda_Be))/CNprop.lambda_Be;
            N26(kk) = N26(kk-1)*exp(-dt*CNprop.lambda_Al) + P26z*(1-exp(-dt*CNprop.lambda_Al))/CNprop.lambda_Al;
    
        end;

        TaBe(i,j) = -1/CNprop.lambda_Be*log(1-N10(ntime)*CNprop.lambda_Be/P10tot);
        
        NBe_M(i,j) = N10(ntime);
        NAl_M(i,j) = N26(ntime);

        
        %figure
        %subplot(3,1,1); hold on; box on;
        %line(time,N10,'color','k'); 
        %title('N10');
        
        %subplot(3,1,2); hold on; box on;
        %line(time,N26,'color','k'); 
        %title('N26');

        %subplot(3,1,3); hold on; box on;
        %line(time,N26./N10,'color','k'); 
        %title('N26/N10');

        
        waitbar(((i-1)*nx+j)/(nx*ny),h);
    
    end;
end;
 
close(h);

 
loaddata;
save TaBe.mat TaBe NBe_M NAl_M bed abrasion noice;

I = find(noice(:) == 1); TaBe(I) = 0;

figure
surf(Xc,Yc,bed,TaBe); shading interp;

figure
surf(Xc,Yc,bed,NAl_M./(NBe_M+1000)); shading interp;

