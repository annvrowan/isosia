function sim2D = karman2d(xs_ext,ys_ext,nu,Lx,seed)

Nx_ext = length(xs_ext);
dx = xs_ext(2)-xs_ext(1);  % dx<>dy does not work yet!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ny_ext = length(ys_ext);
%dy = ys_ext(2)-ys_ext(1);  % dx<>dy does not work yet!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dy=dx;

Nx = 2^nextpow2(Nx_ext);
Nkx = Nx/2 + 1;
dkx = 2*pi/(Nx*dx);
kxs = fftshift(dkx * ((-Nkx+1):(Nkx-2)))';
xs = dx * ((-Nkx+1):(Nkx-2))';
    
Ny = 2^nextpow2(Ny_ext);
Nky = Ny/2 + 1;
dky = 2*pi/(Ny*dy);
kys = fftshift(dky * ((-Nky+1):(Nky-2)))';
ys = dy * ((-Nky+1):(Nky-2))';

kmin = 1/Lx;
    
randn('seed',seed);
beta = 2*nu + 1;
    
P0 = Lx*gamma(nu+0.5)/gamma(nu)/sqrt(pi); % a scaling constant; space domain variance = 1
    
% pre allocations:
P1 = zeros(Nx,Ny);
RHO1  = ones(Nx,Ny)*(1+sqrt(-1));
rho1  = zeros(Nx,Ny);
mask  = zeros(Nx,Ny);
%ks = zeros(Nx,Ny);
    
% the 2-D mask to apply recursively to pick out the band of wavelenghts:
%mask = ( (abs(ks)>=4*dk) & (abs(ks)<=(Nx/4-Nx/32)*dk) );
kx_mask_in = find( abs(fftshift(kxs)) <  4           *dkx );
kx_mask_out= find( abs(fftshift(kxs)) <=  (Nx/4-Nx/32)*dkx );
ky_mask_in = find( abs(fftshift(kys)) <  4           *dky );
ky_mask_out= find( abs(fftshift(kys)) <=  (Ny/4-Ny/32)*dky );
 
kx1 = min(kx_mask_out);
kx2 = min(kx_mask_in);
kx3 = max(kx_mask_in);
kx4 = max(kx_mask_out);
 
ky1 = min(ky_mask_out);
ky2 = min(ky_mask_in);
ky3 = max(ky_mask_in);
ky4 = max(ky_mask_out);
   
%mask = zeros(Nx,Ny);
%mask(kx1:kx4,ky1:ky4) = ones(kx4-kx1+1,ky4-ky1+1);
mask = ones(Nx,Ny);
mask(kx2:kx3,ky2:ky3) = zeros(kx3-kx2+1,ky3-ky2+1);
mask = fftshift(mask);
%image(3*mask);
   
   
% first: save original values dx0 = dx etc.;
dx0  = dx;
xs0  = xs;
dkx0 = dkx;
kxs0 = kxs;
dy0  = dy;
ys0  = ys;
dky0 = dky;
kys0 = kys;
  
% then make the basic series in this sampling: 
for iy=1:Ny
  k_iy2 = (kxs.^2 + kys(iy)^2)/kmin^2;
  ks(:,iy)  = sqrt(k_iy2);
  P1(:,iy) = P0*(1 + k_iy2).^(-(beta+1)/2)*dkx*Nx*dky*Ny;
end
  
% I do not understand completely why "sqrt(Nx)" is required; scaling of fft!
RHO1 = sqrt(Nx*Ny)*sqrt(P1).*(randn(Nx,Ny)+sqrt(-1)*randn(Nx,Ny));
   
rho1 = real(ifft2(RHO1 .* mask));
   
mask = zeros(Nx,Ny);
mask(kx1:kx4,ky1:ky4) = ones(kx4-kx1+1,ky4-ky1+1);
mask(kx2:kx3,ky2:ky3) = zeros(kx3-kx2+1,ky3-ky2+1);
mask = fftshift(mask);
%figure;image(3*mask)
    
for i_cas = 1:2
% next; wavenumbers with dk2 = dk/4;
   dx = dx*Nx/8;
   xs = xs*Nx/8;
    dkx = dkx/(Nx/8);
    kxs = kxs/(Nx/8);
    
    dy = dy*Ny/8;
    ys = ys*Ny/8;
    dky = dky/(Ny/8);
    kys = kys/(Ny/8);
    
    % make series on this coarser grid
    for iy=1:Ny
      k_iy2 = (kxs.^2 + kys(iy)^2)/kmin^2;
      ks(:,iy)  = sqrt(k_iy2);
      P1(:,iy) = P0*(1 + k_iy2).^(-(beta+1)/2)*dkx*Nx*dky*Ny;
    end
    
    % RHO1 is re-used:
    RHO1 = sqrt(Nx*Ny)*sqrt(P1).*(randn(Nx,Ny)+sqrt(-1)*randn(Nx,Ny));
    
    rho_coarse = real(ifft2(RHO1 .* mask));
    
    %Now interpolate:
    iys = find((ys>=min(ys0))&(ys<=max(ys0)));
    iy1 = min(iys)-2; iy2 = max(iys)+2;
    ixs = find((xs>=min(xs0))&(xs<=max(xs0)));
    ix1 = min(ixs)-2; ix2 = max(ixs)+2;
    
    %for iy=iy1:iy2
    %  rho_inty(:,iy-iy1+1) = interp1(xs,rho_coarse(:,iy),xs0);
      rho_inty = interp1(xs(ix1:ix2),rho_coarse(ix1:ix2,iy1:iy2),xs0,'pchip');
    %end
    
    %for ix=1:Nx
      rho_intxy = interp1(ys(iy1:iy2),rho_inty',ys0,'pchip')';
      rho1 = rho1 + rho_intxy;
      
end % loop over cascades in size
%eval(['save c:\tmp\rho',int2str(inu),' rho1 nu']);
% .............................restore original values dx = dx0 etc.;
dx  = dx0;
xs  = xs0;
dkx = dkx0;
kxs = kxs0;
dy  = dy0;
ys  = ys0;
dky = dky0;
kys = kys0;

sim2D = rho1(1:Nx_ext,1:Ny_ext)';
