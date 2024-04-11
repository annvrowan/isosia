function [ice,bed]=show(varargin)

%
% show function for iSOSIA 3
%
% DLE 21/2 2013
%


%load SPM structure
SPM = SPMload;

%default values
fnr = 0;
showtime = 0;
contours = 0;
datatype = ['bed'];
ncon = 50;
withice = 0;
withsedi = 0;
showparticles = 0;
dlim = 0;
tlim = 0;
jpg = 0;
inview = [0,90];
hx = 1;
doclose = 1;
lowdata = -1;
cbar = 1;

%Input variables
if nargin > 0,   
  lArgin = varargin;
  while length(lArgin) >= 1,
    type = lArgin{1};
    lArgin = lArgin(2:end);
    switch type,
     case{'file'}
      fnr = lArgin{1};
      lArgin = lArgin(2:end);
     case{'time'}
      showtime = 1;
     case{'data'}
      datatype = lArgin{1};
      lArgin = lArgin(2:end);
     case{'contours'}
      contours = 1;
     case{'ncon'}
      ncon = lArgin{1};
      lArgin = lArgin(2:end);
     case{'dlim'}
      dlim = 1;
      drange = lArgin{1};
      lArgin = lArgin(2:end);
     case{'tlim'}
      tlim = 1;
      trange = lArgin{1};
      lArgin = lArgin(2:end);
     case{'mindata'}
      lowdata = lArgin{1};
      lArgin = lArgin(2:end);
     case{'view'}
      inview = lArgin{1};
      lArgin = lArgin(2:end);
     case{'hx'}
      hx = lArgin{1};
      lArgin = lArgin(2:end);
     case{'withice'}
      withice = 1;
     case{'withsedi'}
      withsedi = 1;
     case{'particles'}
      showparticles = 1;
     case{'nocbar'}
      cbar = 0;
     case{'noclose'}
      doclose = 0;
     case{'jpg'}
      jpg = 1;
    end;
  end;
end;

%Determine file number
if fnr ~= 0,
  st = load('./status.dat');
  latestfnr = st(3);
  if fnr == 'latest',
    fnr = latestfnr;
    disp(['fnr = ',num2str(fnr)]);
  elseif fnr <= latestfnr,
    fnr = fnr;
  else,
    ['Invalid file number']
    return
  end;
end;

if doclose, 
    close all; 
else 
    figure;
end;
set(gcf,'units','centimeters','position',[10,10,40,25]);
%set(gcf,'color','k','inverthardcopy','off');


%Time series
if showtime,
  
  %load data
  tdata = load('./output/tseries.dat');
  
  subplot(3,3,1);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,3),'-k');
  %plot(tdata(:,1),tdata(:,12),'-r');
  plot(tdata(:,1),tdata(:,25),'-b');
  xlabel('time'); ylabel('H'); title('Mean ice thickness');
  hold off;
  
  subplot(3,3,2);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,4),'-k');
  plot(tdata(:,1),tdata(:,13),'-r');
  plot(tdata(:,1),tdata(:,24),'-b');
  xlabel('time'); ylabel('Nitt'); title('# of velocity itterations');
  hold off;

  subplot(3,3,3);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,2),'-k');
  plot(tdata(:,1),tdata(:,11),'-b');
  xlabel('time'); ylabel('dt'); title('Time step');
  hold off;

  subplot(3,3,4);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,5),'-k');
  xlabel('time'); ylabel('speed'); title('Max average speed');
  hold off;

  subplot(3,3,5);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,6),'-k');
  plot(tdata(:,1),tdata(:,9),'-b');
  %plot(tdata(:,1),tdata(:,15),'-r');
  plot(tdata(:,1),tdata(:,19),'-g');
  plot(tdata(:,1),tdata(:,20),'-y');
  plot(tdata(:,1),tdata(:,21),'-r');
  plot(tdata(:,1),tdata(:,28),'-m');
  xlabel('time'); ylabel('erosion rate'); title('erosion rate');
  hold off;

  
  subplot(3,3,6);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,32)/50,'-k');
  plot(tdata(:,1),tdata(:,33)/50,'-r');
  plot(tdata(:,1),tdata(:,34),'-b');
  plot(tdata(:,1),tdata(:,35)/50,'-g');
  xlabel('time'); ylabel('Right BC'); title('Ice at BC');
  hold off;

  subplot(3,3,7);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,10),'-r');
  xlabel('time'); ylabel('Sealevel T'); title('Sealevel T');
  hold off;

  
  subplot(3,3,8);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,16),'-k');
  plot(tdata(:,1),tdata(:,17),'-b');
  plot(tdata(:,1),tdata(:,18),'-r');
  xlabel('time'); ylabel('Elevation (m)'); 
  title('min, max, and mean elevaion');
  hold off;

  subplot(3,3,9);
  hold on; grid on; box on;
  %plot(tdata(:,1),tdata(:,15),'-k');
  plot(tdata(:,1),tdata(:,29),'-k');
  plot(tdata(:,1),tdata(:,30),'-b');
  plot(tdata(:,1),tdata(:,31),'-r');
  xlabel('time'); ylabel('Particles'); 
  title('Particles');
  hold off;

  
else,
  
  %Retrieve data
  nx = SPM.mesh.nx;
  ny = SPM.mesh.ny;
  L = SPM.mesh.L;
  H = SPM.mesh.H;
  dx = SPM.mesh.dx;
  dy = SPM.mesh.dy;
  include = SPM.data.include;
  fixflag = SPM.data.fixflag;
  
  if fnr == -1,
    
    bed = SPM.data.bed;
    ice = SPM.data.ice;
  
  else,
  
    loaddata;
    %res = sqrt(vxres(:,2:end).^2+vyres(2:end,:).^2);
    %xres = vxres(:,2:end);
    %yres = vyres(2:end,:);
    
  end;
    
  
  %make mesh for plotting
  [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);
  
  %color data
  eval(['data = ',datatype,';']);
  
  %figure
  set(gca,'position',[0.1,0.05,0.7,1],'visible','off');
  hold on; grid on; box on; %axis vis3d;
  load cmap-himal; colormap(map); %colormap
  caxis([0,5]);
  cmap = colormap;
  
  
  if tlim, 
      mintopo = trange(1);
      maxtopo = trange(2);
  else,
      mintopo = min(bed(:));
      maxtopo = max(bed(:));
  end;

  
  
  %if not topography
  if ~strcmp(datatype,'bed'),
    
      %I = find(abs(data(:)) > lowdata);
      if dlim, 
          mindata = drange(1);
          maxdata = drange(2);
          I = find(data(:) > mindata);
      else,
          mindata = min(data(:)); maxdata = max(data(:));
          I = find(abs(data(:)) > 0.000001*(max(data(:))-min(data(:))));
      end; 
      
%       cdata = zeros(size(bed));
%       if (~isempty(I)), cdata(I) = 0.9*(data(I)-mindata)/(maxdata-mindata+1e-16)+4;
%       else cdata(:) = 0.99*(data(:)-mindata)/(maxdata-mindata+1e-6)+4;
%       end; 
      
      cdata =  0.95*(data-mindata)/(maxdata-mindata+1e-16); 
      
      I = find(cdata(:) < 0); cdata(I) = 0;
      I = find(cdata(:) > 1); cdata(I) = 1;
      cdata = cdata + 4; 
      
      zdata = bed+sediment+ice;

      %make Cdata_M
      %Data
      xi = linspace(0,5,length(cmap(:,1)));
      rc = interp1(xi,cmap(:,1),cdata);
      gc = interp1(xi,cmap(:,2),cdata);
      bc = interp1(xi,cmap(:,3),cdata);
      cdata_M(:,:,1) = rc;
      cdata_M(:,:,2) = gc;
      cdata_M(:,:,3) = bc;
  
    

      
  else

        %bed topography
        zdata = bed;
        cdata = 0.99*(zdata-mintopo)/(maxtopo-mintopo);
        I = find(cdata(:) > 0.99); cdata(I) = 0.99;

        %make Cdata_M
        %Elevation
        xi = linspace(0,5,length(cmap(:,1)));
        rc = interp1(xi,cmap(:,1),cdata);
        gc = interp1(xi,cmap(:,2),cdata);
        bc = interp1(xi,cmap(:,3),cdata);
        cdata_M(:,:,1) = rc;
        cdata_M(:,:,2) = gc;
        cdata_M(:,:,3) = bc;
  
        if (withice),
            %Add ice
            ci = 1.1;
            rci = interp1(xi,cmap(:,1),ci);
            gci = interp1(xi,cmap(:,2),ci);
            bci = interp1(xi,cmap(:,3),ci);
            fi = .5*(erf((ice-20)/20)+1);
  
            cdata_M(:,:,1) = (1-fi).*cdata_M(:,:,1)+fi.*rci;
            cdata_M(:,:,2) = (1-fi).*cdata_M(:,:,2)+fi.*gci;
            cdata_M(:,:,3) = (1-fi).*cdata_M(:,:,3)+fi.*bci;
  
            zdata = zdata + ice;       
            
            
            %I = find(ice < 10); zice(I) = bed(I) - 10;
            %cdata = 1.1*ones(size(zice));
            %I = find(ssedi > 0.1); cdata(I) = 2.5;
            %hpi = surf(Xc,Yc,zice,cdata);     
            %set(hpi,'facealpha',.8,'facelighting','none');
            %shading interp;
  
  
        end;
        
        if (withsedi),
                        
            %surface sediment blending
            cs = 2.7;
            rcs = 0.0;
            gcs = 0.0;
            bcs = 0.0;
            
            mins = 0.05;
            ds = 0.05;
    
            fs = .5*(erf((ssedi-mins)/ds)+1);
            I = find(ssedi < 0.01); fs(I) = 0;
  
            cdata_M(:,:,1) = (1-fs).*cdata_M(:,:,1)+fs.*rcs;
            cdata_M(:,:,2) = (1-fs).*cdata_M(:,:,2)+fs.*gcs;
            cdata_M(:,:,3) = (1-fs).*cdata_M(:,:,3)+fs.*bcs;
  
            zsedi = zdata + sediment;

            
            %I = find((sediment < 3)|(ice > 40)); zsedi(I) = bed(I) - 100;
            %cdata = 0.99*(sediment-min(sediment(:)))/(max(sediment(:))-min(sediment(:))+1e-3)+2;
            %hps = surf(Xc,Yc,zsedi,cdata);     
            %set(hps,'facealpha',0.7,'facelighting','none');
            %shading interp;

            
        end;
        
    
        
  
  end;
  
  
        
  %plot surface
  hp = surf(Xc,Yc,zdata,'Cdata',cdata_M);
  set(hp,'facelighting','gouraud','edgelighting','gouraud');
  set(gca,'ambientlightcolor',[.9,.8,.9]);
  %axis equal; 
  set(gca,'dataaspectratio',[1,1,1/hx]);
  shading interp;
  material([.4,.4,.4,2,0.1]);
  view(inview);  
  
  if showparticles,
            
      clear xp;
      clear yp;
      clear zp;
      loadparticles;
      
      
      
      max(zp)
      
      
      %plot3(xp,yp,zp,'.k','markersize',2);
      
  end;
  
  
  
  %keyboard
   
  if contours, 
      
      if tlim, vc = linspace(trange(1),trange(2),ncon);
      elseif withice, vc = linspace(min(bed(:)+sediment(:)+ice(:)),max(bed(:)+sediment(:)+ice(:)),ncon);
      else vc = linspace(min(bed(:)+sediment(:)),max(bed(:)+sediment(:)),ncon);
      end;
  
      if withice,
          contour3(Xc,Yc,bed+sediment+ice+10,vc,'-k'); 
      else, 
          contour3(Xc,Yc,bed+sediment+1,vc,'-k');
      end;
  end;
    
  if cbar, 
      vt = [0,0.25,0.5,0.75,1.0];
      cb = colorbar('position',[0.85,0.1,0.02,0.4]);
      set(cb,'ylim',[0,1],'ytick',vt,'yticklabel',mintopo+vt*(maxtopo-mintopo));
      set(get(cb,'ylabel'),'string',['elevation (m)']);
  end;
 
  if ~strcmp(datatype,'bed'),
    vt = [0,0.25,0.5,0.75,1.0];
    cb = colorbar('position',[0.85,0.55,0.02,0.4]);
    set(cb,'ylim',[4,5],'ytick',vt+4,'yticklabel',mindata+vt*(maxdata-mindata));
    set(get(cb,'ylabel'),'string',datatype);
  end;
  
  light('position',[100,0,500]*1e3);
  %lightangle(inview(1)+0,65); light;
  
  
 
  
  hold off;
  
  if jpg,
    pause(.2);
    
    eval(['print -djpeg90 -r100 ./flic/flic',sprintf('%04d',fnr),'.jpg']);
    
  end;

  if SPM.mesh.docelldata,
    load input/nodelist.mat;
    hold on;
    for i=1:length(nnx),
      plot3(Xc(nny(i),nnx(i)),Yc(nny(i),nnx(i)),3000* ...
            ones(size(nnx)),'.r','markersize',30);
      end;
  end;
  
end;

