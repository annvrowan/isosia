function burial_path(fnr)

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
Pspal_Be = 3.9; %atoms/(g*yr)
Pnmc_Be = 0.0969;
Pfm_Be = 0.0851;

Pspal_Al = 26.8; %atoms/(g*yr)
Pnmc_Al = 0.820;
Pfm_Al = 0.703;

Lspal_Be = 160; %g/cm^2
Lspal_Al = 160; %g/cm^2
TBe = 1.39e6;
TAl = 0.705e6;
lambda_Be = log(2)/TBe;
lambda_Al = log(2)/TAl;
rho = 2.65; %density g/cm3
            
%************* dO18 *****************
load lisiecki.mat
I = find(Age < 2.7);
dO = d18O(I);
time = Age(I);
ntime = length(time);
dt = diff(time);

dfis = 4.2;
I = find(dO > dfis);
Tfis = sum(dt(I))
dOfis = dfis*ones(size(dO));
Is = find(dO > dfis); dOfis(Is) = dO(Is);
itime = linspace(Tfis,0,Nt);
fitime = linspace(Tfis,0,length(Is));

%hold on; box on;
%set(gca,'position',[.05,.1,.9,.85]);
%patch(time,dOfis,[63,160,203]/256);
%line(time,dO,'color','k','linewidth',.2);
%set(gca,'ydir','reverse');
%set(gca,'xtick',[0:0.5:2.5]);
%set(gca,'ytick',[3:1:5]);
%axis([0,2.7,3,5.2]);
%**************************************

TaBe = zeros(ny,nx);
NBe_M = zeros(ny,nx);
NAl_M = zeros(ny,nx);
h = waitbar(0,'Computing...');
for i=1:ny,
    for j=1:nx,
 
        
        ipath = zeros(1,Nt);
        ipath(:) = burial(i,j,:);
        fipath = interp1(itime,ipath,fitime);

        fpath = zeros(size(time));
        pfac = ones(size(time));
        fpath(Is) = fliplr(fipath);
        pfac(Is) = 0;
        mb = 0;
        for kk=1:length(time),
            if fpath(kk) == 0,
                fpath(kk) = mb;
            elseif fpath(kk) > mb,
                mb = fpath(kk);
            end;
        end;

        %allocate
        NBe = zeros(ntime,1);
        NAl = zeros(ntime,1);

        for kk=(ntime-1):-1:1,
    
            dt = (time(kk+1)-time(kk))*1e6;
            Pz_Be = pfac(kk)*(Pspal_Be+Pnmc_Be+Pfm_Be)*exp(-rho*100*fpath(kk)/Lspal_Be);
            Pz_Al = pfac(kk)*(Pspal_Al+Pnmc_Al+Pfm_Al)*exp(-rho*100*fpath(kk)/Lspal_Al);
            NBe(kk) = NBe(kk+1)*exp(-dt*lambda_Be) + Pz_Be*(1-exp(-dt*lambda_Be))/lambda_Be;
            NAl(kk) = NAl(kk+1)*exp(-dt*lambda_Al) + Pz_Al*(1-exp(-dt*lambda_Al))/lambda_Al;
    
        end;

        TaBe(i,j) = -1/lambda_Be*log(1-NBe(1)*lambda_Be/(Pspal_Be+ ...
                                                       Pnmc_Be+Pfm_Be));
        NBe_M(i,j) = NBe(1);
        NAl_M(i,j) = NAl(1);

        waitbar(((i-1)*nx+j)/(nx*ny),h);
    
    end;
end;
 
close(h);

%figure(2);
%hold on; grid on; box on;
%line(itime,ipath,'color','r');
%line(fitime,fipath,'color','b');
%line(time,fpath,'color','k');

%figure(3);
%hold on; grid on; box on;
%line(time,NBe,'color','r');
%line(time,NAl,'color','b');
 
loaddata;
save TaBe.mat TaBe NBe_M NAl_M bed abrasion;



surf(Xc,Yc,bed,TaBe); shading interp;

figure
surf(Xc,Yc,bed,NAl_M./NBe_M); shading interp;

