function flowlines()

close all;


%load input
load data_for_flowlines.mat
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;

%step size along flowline
dd = 0.25*(dx+dy);

%minimum ice velocity
minvelo = 1e-8;

%cosmo parameters
Pspal_Be = 30; %atoms/(g*yr)
Lspal = 150; %g/cm^2
TBe = 1.378e6;
lambda_Be = log(2)/TBe;
rhoi = 0.9; %density g/cm3
NBe0 = 0; %starting Be concentration 

%top surface
ti = bed + ice + sediment;

%make mesh
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

%find margin cells
[mrow,mcol] = find(icemargin == 1);
Nm = length(mrow);

%loop margin cells
for i=1:Nm,
    
    %initialize
    lx = [];
    ly = [];
    lb = [];
    lz = [];
    lbf = [];
    ltime = [];
    ldist = [];
    lNBe = [];
    time = 0;
    dist = 0;
    i
    
    %margin cell is starting point
    nn = 1; hit = 1;
    ni = mrow(i);
    nj = mcol(i);
    xx = Xc(ni,nj);
    yy = Yc(ni,nj);
    bz = 0.0;
    bf = 0.0; %surface
    NBe = NBe0;
    
    %save starting position to flowline
    lx(nn) = xx;
    ly(nn) = yy;
    lb(nn) = bz;  %burial depth
    lbf(nn) = bf; %burial ratio
    lz(nn) = ti(ni,nj)-bz;
    ltime(nn) = time;
    ldist(nn) = dist;
    lNBe(nn) = NBe;
    
    %cell midpoint velocities
    vxx = .5*(vx(ni,nj) + vx(ni,nj+1));
    vyy = .5*(vy(ni,nj) + vy(ni+1,nj));
    vv = sqrt(vxx^2 + vyy^2);
        
    %project next position along flowline
    if (vv > minvelo),
        dt = dd/vv;
        ddx = vxx*dt;
        ddy = vyy*dt;
        xx = xx - ddx;
        yy = yy - ddy;
        time = time + dt;
        dist = dist + dd;
        Pz = Pspal_Be*exp(-rhoi*100*bz/Lspal);
        NBe = NBe*exp(-dt*lambda_Be) + Pz*(1-exp(-dt*lambda_Be))/lambda_Be;
    else
        hit = -1;
    end;
     
    itt = 0;
    %loop until no ice or outside grid
    while (hit > 0)&(itt < 4000),%nx+ny),    
    
        %find cell row and col numbers of new point
        ni = ceil(yy/dy);
        nj = ceil(xx/dx);
        
        %if outside grid
        if ((ni > (ny-1))|(nj > (nx-1))|(ni < 1)|(nj < 1)),
            %if (ni >= ny), ni = ny-1; end;
            %if (nj >= nx), nj = nx-1; end;
            hit = -1;
            %elseif icemargin(ni,nj) == 1,
            %hit = -1;
        else,
            
            %increase index
            nn = nn + 1;

            %update burial
            bz = bf*ice(ni,nj);
            vzz = (1-bf)*massb(ni,nj)-bf*bmelt(ni,nj);
            ddz = vzz*dt;
            bz = bz + ddz;
            if (bz < 0), 
                bz = 0; 
            elseif (bz > ice(ni,nj)),
                bz = ice(ni,nj);
            end;
            if (ice(ni,nj) > 1),
                bf = bz/ice(ni,nj);
            else bf = 0;
            end;
            
            %get horizontal velocity components by bilinear interpolation
            yyn = (yy - (ni-1)*dy)/dy;
            xxn = (xx - (nj-1)*dx)/dx;
            vxx = [1-xxn,xxn]*[vx(ni,nj),vx(ni+1,nj);vx(ni,nj+1),vx(ni+1,nj+1)]*[1-yyn;yyn];
            vyy = [1-xxn,xxn]*[vy(ni,nj),vy(ni+1,nj);vy(ni,nj+1),vy(ni+1,nj+1)]*[1-yyn;yyn];
                        
            %ice surface
            hi = [1-xxn,xxn]*[ti(ni,nj),ti(ni+1,nj);ti(ni,nj+1),ti(ni+1,nj+1)]*[1-yyn;yyn];
            
            %save positions
            lx(nn) = xx;
            ly(nn) = yy;
            lb(nn) = bz;
            lbf(nn) = bf;
            lz(nn) = hi-bz;
            ltime(nn) = time;
            ldist(nn) = dist;
            lNBe(nn) = NBe;
           
            %ice speed
            vv = sqrt(vxx^2 + vyy^2);
        
            %project next step
            if (vv > minvelo),
                dt = dd/vv;
                ddx = vxx*dt;
                ddy = vyy*dt;
                ddz = vzz*dt;
                xx = xx - ddx;
                yy = yy - ddy;
                time = time + dt;
                dist = dist + dd;
                Pz = Pspal_Be*exp(-rhoi*100*bz/Lspal);
                NBe = NBe*exp(-dt*lambda_Be) + Pz*(1-exp(-dt*lambda_Be))/lambda_Be;
            else
                hit = -1;
            end;
            
    
           
        end;
        itt = itt + 1;
    end;
    
    
    flines{i}.x = lx;
    flines{i}.y = ly;
    flines{i}.z = lz;
    flines{i}.bz = lb;
    flines{i}.bf = lbf;
    flines{i}.time = ltime;
    flines{i}.dist = ldist;
    flines{i}.NBe = lNBe;


    plot3(lx,ly,lz)    
    hold on;    
    
end;

save flines.mat flines


