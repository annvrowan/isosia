function [bed,ice,sediment,Vs]=show(varargin)

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
dlim = 0;
jpg = 0;
inview = [0,90];
hx = 1;
doclose = 1;

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
     case{'view'}
      inview = lArgin{1};
      lArgin = lArgin(2:end);
     case{'hx'}
      hx = lArgin{1};
      lArgin = lArgin(2:end);
     case{'withice'}
      withice = 1;
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
  plot(tdata(:,1),tdata(:,12),'-r');
  %plot(tdata(:,1),tdata(:,25),'-b');
  xlabel('time'); ylabel('H'); title('Mean ice thickness');
  hold off;
  
  subplot(3,3,2);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,4),'-k');
  plot(tdata(:,1),tdata(:,13),'-r');
  %plot(tdata(:,1),tdata(:,24),'-b');
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
  plot(tdata(:,1),tdata(:,15),'-r');
  plot(tdata(:,1),tdata(:,19),'-g');
  plot(tdata(:,1),tdata(:,20),'-m');
  plot(tdata(:,1),tdata(:,21),'-y');
  xlabel('time'); ylabel('erosion rate'); title('erosion rate');
  hold off;

  
  subplot(3,3,6);
  hold on; grid on; box on;
  plot(tdata(:,1),tdata(:,7),'-k');
  plot(tdata(:,1),tdata(:,8),'-r');
  plot(tdata(:,1),tdata(:,14),'-b');
  xlabel('time'); ylabel('sediment transfer'); title('Sediment transport');
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
  plot(tdata(:,1),tdata(:,23),'-k');
  xlabel('time'); ylabel('Sediment (m)'); 
  title('mean sediment');
  hold off;

  
else,
  
  %Retrieve data
  nx = SPM.mesh.nx;
  ny = SPM.mesh.ny;
  L = SPM.mesh.L;
  H = SPM.mesh.H;
  dx = SPM.mesh.dx;
  dy = SPM.mesh.dy;
  maxacc = SPM.data.maxacc;
  fixflag = SPM.data.fixflag;
  
  if fnr == -1,
    
    bed = SPM.data.bed;
    ice = SPM.data.ice;
  
  else,
  
    loaddata;
        
  end;
    
  avvelo = mean(abs(sliding(:)+deformation(:)))
  avero = mean(abrasion(:)+quarrying(:)+frost(:)+hillslope(:)+landslide(:))
  
  %make mesh for plotting
  [Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

  %bed topography information
  minbed = 100*floor(min(bed(:))/100);
  maxbed = 100*ceil(max(bed(:))/100);
  
  
  %color data
  eval(['data = ',datatype,';']);
  
  
  %figure
  set(gca,'position',[0.1,0.05,0.7,1],'visible','off');
  hold on; grid on; box on; %axis vis3d;
  load cmap-himal; colormap(map); %colormap
  
  %bed topography
  cdata = 0.99*(bed-min(bed(:)))/(max(bed(:))-min(bed(:)));
  zdata = bed+sediment;
  cdata = 0.99*(zdata-min(zdata(:)))/(max(zdata(:))-min(zdata(:)));
    
  %if only topography
  if strcmp(datatype,'bed'),

    if withice,
      I = find(ice > 10);
      cdata(I) = 0.99*(ice(I)-min(ice(I)))/(max(ice(I))-min(ice(I))+1e-6)+1;
      zdata = bed+ice;
    
    end
    
  else, %if (strcmp(datatype,'ice')|strcmp(datatype,'te2')|strcmp(datatype,'sliding')|strcmp(datatype,'deformation')|strcmp(datatype,'quarrying')|strcmp(datatype,'pw2')|strcmp(datatype,'bmelt')),
    
    %I = find(abs(data(:)) > 0.0001*(max(data(:))-min(data(:))));
    I = find(abs(data(:)) > -1);
    if dlim, 
      mindata = drange(1);
      maxdata = drange(2);
    else,
      mindata = min(data(:)); maxdata = max(data(:));
      %mindata = min(data(I)); maxdata = max(data(I));
      %order = floor(log10(max(abs(maxdata)+eps,abs(mindata)+eps)));
      %mindata = 10^order*floor(mindata)/10^order;
      %maxdata = 10^order*ceil(maxdata)/10^order;
    end;
    if (~isempty(I)), cdata(I) = 0.99*(data(I)-mindata)/(maxdata-mindata+1e-16)+4;
    else cdata(:) = 0.99*(data(:)-mindata)/(maxdata-mindata+1e-6)+ ...
          4;
    end;
    zdata = bed+ice+sediment;
    
      
  %else,

  %  if dlim, 
  %    mindata = drange(1);
  %    maxdata = drange(2);
  %  else,
  %    mindata = min(data(:)); maxdata = max(data(:));
  %    order = floor(log10(max(abs(maxdata),abs(mindata))));
  %    mindata = 10^order*floor(mindata)/10^order;
  %    maxdata = 10^order*ceil(maxdata)/10^order;
  %  end;
  %  cdata(:) = 0.99*(data(:)-mindata)/(maxdata-mindata+1e-6)+4;
  %  zdata = bed+ice;
        
  end;
      
  %plot bed topography
  hp = surf(Xc,Yc,zdata,cdata);
  set(hp,'facelighting','gouraud','edgelighting','gouraud');
  set(gca,'ambientlightcolor',[.9,.8,.9]);
  %axis equal; 
  set(gca,'dataaspectratio',[1,1,1/hx]);
  shading interp;
  caxis([0,5])
  material([.4,.4,.4,2,0.1]);
  %view(inview);  
  
  %keyboard
  
  if contours, 
    if withice,
      contour3(Xc,Yc,bed+sediment+ice+10,ncon,'-k'); 
    else,
      contour3(Xc,Yc,bed+sediment+1,ncon,'-k');
    end;
  end;
  
  
  vt = [0,0.25,0.5,0.75,1.0];
  cb = colorbar('position',[0.85,0.1,0.02,0.4]);
  set(cb,'ylim',[0,1],'ytick',vt,'yticklabel',minbed+vt*(maxbed-minbed));
  set(get(cb,'ylabel'),'string',['elevation (m)']);
  
 
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

% debris = Vs{1};
% debris(debris<0.01) = NaN;
% mean_debris = mean(mean(debris,'omitnan'),'omitnan');
% std_debris = std(std(debris,'omitnan'),'omitnan');


% ela = mean(bed(meltrate==0 & ice > 0));
% ela;

dh = Vs{1};
dh(dh>5)=5;
dh = flipud(dh);
velo = sliding+deformation;
velo(velo>200)=200;
velo = flipud(velo);
ice(ice<0)=0;
ice = flipud(ice);

%save('ice','ice','-ascii')
%save('dh','dh','-ascii')
save('velo','velo','-ascii')




% 
% meltrate(ice==0)=NaN;
% sum_meltrate = sum(sum(meltrate,'omitnan'),'omitnan')
% std_meltrate = std(std(meltrate,'omitnan'),'omitnan')






