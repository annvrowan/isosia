function movieframe(fnrs)

close all
set(gcf,'units','pixels','position',[100,100,1920,1080]);
set(gcf,'color',.0*[1,1,1],'paperpositionmode','auto','inverthardcopy','off');
set(gca,'position',[0,0,1,1],'visible','off');
axis([0,1,0,1]); hold on;


%user parameters
mintopo = -2000;
maxtopo = 4000;
hx = 20;
inview = [-120,16];
%inview = [0,90];
Dc = 2500;
maxfnr = 740;
showsliding = 0;
showqw = 0;
showte = 0;
showice = 1;
showabrasion = 0;
showelevation = 1;
showTaBe = 0;
showslidinglines = 0;
timeseries = 0;
flowlines = 0;

minsliding = 0;
maxsliding = 10;
minqw = 0;
maxqw = 0.01;
minte = 0;
maxte = 200;
minabrasion = 0;
maxabrasion = 400;
minTaBe = 20;
maxTaBe = 160;

%load data
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
L = SPM.mesh.L;
H = SPM.mesh.H;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;

%make mesh for plotting
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

%colormap
load cmap; colormap(map); %colormap
caxis([0,5]);
cmap = colormap;
col1 = [215,211,167]/256;
col2 = [100,157,183]/256;

%curve legend
%line([0.7,0.75],[0.875,0.875],'color',col1,'linewidth',3);
%line([0.7,0.75],[0.825,0.825],'color',col2,'linewidth',3);
%plot(0.75,0.825,'o','markersize',15,'markerfacecolor',col2,'markeredgecolor',col2);
%text(0.76,0.875,'Sea level Temperature','HorizontalAlignment','left','VerticalAlignment','middle','color',col1,'fontsize',32);
%text(0.76,0.825,'Ice volume','HorizontalAlignment','left','VerticalAlignment','middle','color',col2,'fontsize',32);

ax1 = axes('position',[0.01,-0.075,.98,1]);
set(gca,'xlim',[0,L],'ylim',[0,H],'zlim',[-4000,6000]);
set(gca,'dataaspectratio',[1,1,1/hx]);
set(gca,'visible','off');
view(inview); camproj(gca,'perspective');  
hold on;
ax1 = gca;

%make timeseries panel
if timeseries,
    pos = get(gcf,'position');
    ax2 = axes('position',[0.05,0.675,0.55,0.2],'visible','on');
    tsc = 10*300/265; %time scaling
    hold on; box on;
    set(gca,'color','none');
    axis([0,300e3,-1,7]);
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    set(gca,'xcolor',col1,'ycolor',col1);
    set(gca,'fontsize',30,'ytick',[0,3,6],'xtick',[0:100e3:400e3],'xticklabel',[0:100:300]);
    xlabel('Time (kyr)'); ylabel('Temperature');
end;

cb = colorbar('horizontal','position',[.7,.675,.25,.03]);
set(cb,'color',col1,'LimitsMode','manual','TickLabelsMode','manual','TicksMode','manual');



if showsliding,
    
    set(cb,'limits',[0.8,1.0],'ticks',[0.8:0.05:1.0],'ticklabels', ...
           [0:10:40],'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Sliding speed (m/y)','fontsize',32,'color',col1);

elseif showqw,
    
    nt = 3;
    set(cb,'limits',[0.8,1.0],'ticks',linspace(0.8,1.0,nt),'ticklabels', ...
           linspace(minqw,maxqw,nt),'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Water flux (m2/y)','fontsize',32,'color',col1);

elseif showte,
    
    nt = 3;
    set(cb,'limits',[0.8,1.0],'ticks',linspace(0.8,1.0,nt),'ticklabels', ...
           linspace(minte,maxte,nt),'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Normalized pressure (m)','fontsize',32,'color',col1);

elseif showelevation,
    
    set(cb,'limits',[0.0,0.2],'ticks',[0.0:0.2/3:0.2],'ticklabels', ...
           [-2000:2000:4000],'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Elevation (m a.s.l.)','fontsize',32,'color',col1);


elseif showabrasion
    
    set(cb,'limits',[0.8,1.0],'ticks',[0.8:0.05:1.0],'ticklabels', ...
           [0:100:400],'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Erosion (m)','fontsize',32,'color',col1);
    
elseif showTaBe,
    
    set(cb,'limits',[0.8,1.0],'ticks',[0.8:0.05:1.0],'ticklabels', ...
           [0:40:160],'color',col1);
    set(cb,'fontsize',30);
    xlabel(cb,'Apparent exposure age (kyr)','fontsize',32,'color',col1);
    

end;

%time series info
if timeseries,
    time = SPM.mprop.Temp(1,:);
    Temp = SPM.mprop.Temp(2,:);
    
    time = linspace(0,100e3,1000);
    Temp = 2*ones(size(time));
    
    %load time series data
    tdata = load('./output/tseries.dat');
    [ttime,Iu] = unique(tdata(:,1));
    Icevol = interp1(tdata(Iu,1),tdata(Iu,3),time);
end

  
for i=1:length(fnrs),
    
    fnr = fnrs(i);
    loaddata;

    
    if showelevation,
    
        %********* make bed surface ************* 
        %Elevation
        zdata = bed;
        cdata = 0.99*(zdata-mintopo)/(maxtopo-mintopo);
        I = find(cdata(:) > 0.99); cdata(I) = 0.99;
        I = find(cdata(:) < 0.01); cdata(I) = 0.01;

    elseif showabrasion,
        
        zdata = bed;
        cdata =  0.95*(abrasion-minabrasion)/(maxabrasion-minabrasion);      
        I = find(cdata(:) < 0); cdata(I) = 0;
        I = find(cdata(:) > 1); cdata(I) = 1;
        cdata = cdata + 4; 
        
    elseif showsliding, 
    
        zdata = bed;
        cdata =  0.95*(sliding-minsliding)/(maxsliding-minsliding);      
        I = find(cdata(:) < 0); cdata(I) = 0;
        I = find(cdata(:) > 1); cdata(I) = 1;
        cdata = cdata + 4; 
            
    elseif showqw, 
    
        zdata = bed;
        cdata =  0.95*(qw-minqw)/(maxqw-minqw);      
        I = find(cdata(:) < 0); cdata(I) = 0;
        I = find(cdata(:) > 1); cdata(I) = 1;
        cdata = cdata + 4; 
    
    elseif showte, 
    
        zdata = bed;
        cdata =  0.95*(te-minte)/(maxqw-minte);      
        I = find(cdata(:) < 0); cdata(I) = 0;
        I = find(cdata(:) > 1); cdata(I) = 1;
        cdata = cdata + 4; 
    
    elseif showTaBe, 
    
        %load exposure ages
        load TaBe.mat;
        TaBe = TaBe*1e-3;
        
        zdata = bed;
        cdata =  0.95*(TaBe-minTaBe)/(maxTaBe-minTaBe);      
        I = find(cdata(:) < 0); cdata(I) = 0;
        I = find(cdata(:) > 1); cdata(I) = 1;
        cdata = cdata + 4; 
    
    end;

    %make colormap
    xi = linspace(0,5,length(cmap(:,1)));
    rc = interp1(xi,cmap(:,1),cdata);
    gc = interp1(xi,cmap(:,2),cdata);
    bc = interp1(xi,cmap(:,3),cdata);
    cdata_M(:,:,1) = rc;
    cdata_M(:,:,2) = gc;
    cdata_M(:,:,3) = bc;
    
    %plot bed surface
    axes(ax1); cla;
    hp = surf(Xc,Yc,zdata,'Cdata',cdata_M); caxis('manual');
    set(hp,'facelighting','gouraud','edgelighting','gouraud');
    set(gca,'ambientlightcolor',[.9,.8,.9]);
    light
    shading interp;
    material([.4,.4,.4,2,0.1]);
    
    %ice color
    ci = 1.1;
    rci = interp1(xi,cmap(:,1),ci);
    gci = interp1(xi,cmap(:,2),ci);
    bci = interp1(xi,cmap(:,3),ci);
    
    
    if showice,
        
        fi = .5*(erf((ice-100)/20)+1);
    
        cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*rci;
        cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*gci;
        cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*bci;
       
        if showsliding == 1,
        
            %make data color
            cdata =  0.95*(sliding-minsliding)/(maxsliding-minsliding);      
            I = find(cdata(:) < 0); cdata(I) = 0;
            I = find(cdata(:) > 1); cdata(I) = 1;
            cdata = cdata + 4; 
    
            xi = linspace(0,5,length(cmap(:,1)));
            rc = interp1(xi,cmap(:,1),cdata);
            gc = interp1(xi,cmap(:,2),cdata);
            bc = interp1(xi,cmap(:,3),cdata);
            cdata_S(:,:,1) = rc;
            cdata_S(:,:,2) = gc;
            cdata_S(:,:,3) = bc;
    
            fi = .5*(erf((sliding-.5)/.025)+1);

            fi = ones(size(sliding));
            
            %cdata_M = (1-fi)*cdata_M + fi*cdata_S;
            cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*cdata_S(:,:,1);
            cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*cdata_S(:,:,2);
            cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*cdata_S(:,:,3);
       
        elseif showqw == 1,
        
            %make data color
            cdata =  0.95*(qw-minqw)/(maxqw-minqw);      
            I = find(cdata(:) < 0); cdata(I) = 0;
            I = find(cdata(:) > 1); cdata(I) = 1;
            cdata = cdata + 4; 
    
            xi = linspace(0,5,length(cmap(:,1)));
            rc = interp1(xi,cmap(:,1),cdata);
            gc = interp1(xi,cmap(:,2),cdata);
            bc = interp1(xi,cmap(:,3),cdata);
            cdata_S(:,:,1) = rc;
            cdata_S(:,:,2) = gc;
            cdata_S(:,:,3) = bc;
    
            fi = .5*(erf((qw-.001)/.001)+1);
         
            %cdata_M = (1-fi)*cdata_M + fi*cdata_S;
            cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*cdata_S(:,:,1);
            cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*cdata_S(:,:,2);
            cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*cdata_S(:,:,3);
       
        elseif showte == 1,
        
            %make data color
            cdata =  0.95*(te-minte)/(maxte-minte);      
            I = find(cdata(:) < 0); cdata(I) = 0;
            I = find(cdata(:) > 1); cdata(I) = 1;
            cdata = cdata + 4; 
    
            xi = linspace(0,5,length(cmap(:,1)));
            rc = interp1(xi,cmap(:,1),cdata);
            gc = interp1(xi,cmap(:,2),cdata);
            bc = interp1(xi,cmap(:,3),cdata);
            cdata_S(:,:,1) = rc;
            cdata_S(:,:,2) = gc;
            cdata_S(:,:,3) = bc;
    
            fi = .5*(erf((te-10)/10)+1);
         
            %cdata_M = (1-fi)*cdata_M + fi*cdata_S;
            cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*cdata_S(:,:,1);
            cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*cdata_S(:,:,2);
            cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*cdata_S(:,:,3);
       
          
        end;
            
        zdata = zdata + ice - 100;       
        hp = surf(Xc,Yc,zdata,'Cdata',cdata_M); caxis('manual');
        set(hp,'facelighting','gouraud','edgelighting','gouraud'); 
        %light
        set(hp,'facealpha',.8);
        %set(hp,'facealpha',.5);
        shading interp;
        material([.6,.3,.4,2,1]);
        
        if showslidinglines,
        
            lcolor = .4*[1,1,1];
            lw0 = 1;
            msliding = mean(sliding,2);
            mx = floor(.9*nx);
            for i=10:10:(ny-10),
                dxs = 5e3*msliding(i);
                enx = mx-round(dxs/dx);
                lw = lw0 + 2*round(dxs/dx)/nx;
                if (enx < 1), enx = 1; end;
                line([Xc(i,mx),Xc(i,mx)-dxs],[Yc(i,mx),Yc(i,mx)], ...
                     [zdata(i,mx),zdata(i,enx)],'color',lcolor,'linewidth',lw);
            end;
            
        end;
        

    end

    %******* patch sides ***********

    bc = [.85,.85,.8];
    
    yp = [Yc(:,1)',max(Yc(:)),min(Yc(:))];
    zp = [bed(:,1)',-Dc,-Dc];
    patch(min(Xc(:))*ones(size(yp)),yp,zp,bc);

    xp = [Xc(1,:),max(Xc(:)),min(Xc(:))];
    zp = [bed(end,:),-Dc,-Dc];
    patch(xp,max(Yc(:))*ones(size(xp)),zp,bc);

    
    
    %time series
    if timeseries,
        mtime = SPM.mesh.filetime*fnr;
   
        axes(ax2); cla;
        line(time*tsc,Temp,'color',col1,'linewidth',3);
    
        %I = find(time <= mtime);
        %line(time(I)*tsc,Icevol(I)/150,'color',col2,'linewidth',3);
        %dp = interp1(time,Icevol,mtime);
        %plot(mtime*tsc,dp/150,'o','markersize',15,'markerfacecolor',col2,'markeredgecolor',col2);
        
        dp = interp1(time,Temp,mtime);
        plot(mtime*tsc,dp,'o','markersize',15,'markerfacecolor',col1,'markeredgecolor',col1);
    end;

    
    if flowlines,
        addflowlines;
    end;
    
        
        
    pause(0.25);
    
    
    eval(['print -djpeg90 -r0 ./flic/flic',num2str(fnr),'.jpg']);
    
end;

