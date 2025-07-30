function no = noise(ix,iy,my,L,X,Y,feed)
%function no = noise(ix,iy,my,L,x,y,feed)

no = zeros(size(X));

%von Karman noise on regular grid
noise = karman2d(ix(:),iy(:),my,L,feed);

%interpolate onto mesh
%inoise = interp2(ix,iy,noise,x,y,'nearest');
%inoise = interp2(ix,iy,noise,x,y);

no(:) = noise;

%normalize noise
no = no/(max(no(:))-min(no(:))) - min(no(:));
