function hypsometry(varargin)


% hypsometry - shows hypsometry for SPM model
%
%   DLE 11/3 2013
%

close all;

set(gcf,'units','centimeters','position',[10,10,20,25]);
set(gcf,'color','w');

%load SPM structure
SPM = SPMload;
nx = SPM.mesh.nx;
ny = SPM.mesh.ny;
Nc = nx*ny;

%default values
fnr = 0;
Nh = 30;


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

%initial bed elevation
bed0 = SPM.data.bed;

if fnr == 0,
  
  bed = bed0;
  
else,

  %load data
  fid = fopen(['./output/output',num2str(fnr),'.dat']);
  xc = fread(fid,[ny,nx],'double');
  yc = fread(fid,[ny,nx],'double');
  bed = fread(fid,[ny,nx],'double');
  fclose(fid);

end;

maxele = max([bed(:);bed0(:)]);
minele = min([bed(:);bed0(:)]);

edges = linspace(minele,maxele,Nh+1);
N0 = histc(bed0(:),edges)/Nc;
N = histc(bed(:),edges)/Nc;

hold on; box on;
for i=1:Nh,
  patch([0,N(i),N(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],[.3,.4,.6]);
end;
%for i=1:Nh,
%  line([0,N0(i),N0(i),0],[edges(i),edges(i),edges(i+1),edges(i+1)],'color','b');
%end;

set(gca,'xlim',[0,.15],'ylim',[-1000,2500]);


print -dtiff -r100 hypso60.tif