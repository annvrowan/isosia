function [ice,hw,bed]=makeframe(varargin)

%
% show function for iSOSIA 3 updated to function called makeframe
% to visualise results from later versions of spm (spm-3.4.0 onwards)
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
inview = [10,45];
hx = 1;
flowlines = 0;
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
     case{'flowlines'}
      flowlines = 1;
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


%Retrieve data
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
L = SPM.mesh.L;
H = SPM.mesh.H;
dx = SPM.mesh.dx;
dy = SPM.mesh.dy;
fixflag = SPM.data.fixflag;
  
if fnr == -1,
    
    bed = SPM.data.bed;
    ice = SPM.data.ice;
    
else,
  
    loaddata;
        
end;
    
%make mesh for plotting
[Xc,Yc] = meshgrid([1:nx]*dx-dx/2,[1:ny]*dy-dy/2);

%bed topography information
%minbed = 100*floor(min(bed(:))/100);
%maxbed = 100*ceil(max(bed(:))/100);
  
%color data
eval(['data = ',datatype,';']);
  
%figure
set(gca,'position',[0.1,0.05,0.7,1],'visible','off');
hold on; grid on; box on; %axis vis3d;
load cmap-himal; colormap(map); %colormap
  
%bed topography
cdata = 0.99*(bed-min(bed(:)))/(max(bed(:))-min(bed(:)));
zdata = bed+sediment;
if dlim, 
    mindata = drange(1);
    maxdata = drange(2);
else,
    mindata = min(zdata(:));
    maxdata = max(zdata(:));
end;
  
%cdata = 0.99*(zdata-mindata)/(maxdata-mindata);
    
%if only topography
if strcmp(datatype,'bed'),

    if withice,
        I = find(ice > 10);
        cdata(I) = 0.99*(ice(I)-min(ice(I)))/(max(ice(I))-min(ice(I))+1e-6)+1;
        zdata = bed+ice;
    
    end
    
else, %if (strcmp(datatype,'ice')|strcmp(datatype,'te2')|strcmp(datatype,'sliding')|strcmp(datatype,'deformation')|strcmp(datatype,'quarrying')|strcmp(datatype,'pw2')|strcmp(datatype,'bmelt')),
    
    %I = find(abs(data(:)) > 0.000001*(max(data(:))-min(data(:))));
    I = find(ice(:) > 10);
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
    else cdata(:) = 0.99*(data(:)-mindata)/(maxdata-mindata+1e-6)+4;
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
view(inview);  

%keyboard
   
if dlim, vc = linspace(drange(1),drange(2),ncon);
elseif withice, vc = linspace(min(bed(:)+sediment(:)+ice(:)),max(bed(:)+sediment(:)+ice(:)),ncon);
else vc = linspace(min(bed(:)+sediment(:)),max(bed(:)+sediment(:)),ncon);
end;

if contours, 
    if withice,
        contour3(Xc,Yc,bed+sediment+ice+10,vc,'-k'); 
    else, 
        contour3(Xc,Yc,bed+sediment+1,vc,'-k');
    end;
end;

if flowlines,
    addflowlines;
end;


  
%vt = [0,0.25,0.5,0.75,1.0];
%cb = colorbar('position',[0.85,0.1,0.02,0.4]);
%set(cb,'ylim',[0,1],'ytick',vt,'yticklabel',mindata+vt*(maxdata-mindata));
%set(get(cb,'ylabel'),'string',['elevation (m)']);


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
  

