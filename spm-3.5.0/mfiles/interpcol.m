function [col] = interpcol(minval,maxval,val)

% Function to interpolate the jet colormap given a data value range
% (min, max) and a data point vector (val)


% Prepare colormap
y = [0:0.0625:1]';
ry = flipud(y);

cmap = zeros(64,3);

cmap(24:40,1) = y;
cmap(41:55,1) = 1;
cmap(56:64,1) = ry(1:9);

cmap(8:24,2) = y;
cmap(24:39,2) = 1;
cmap(40:56,2) = ry;

cmap(1:8,3) = y(10:17);
cmap(9:23,3) = 1;
cmap(24:40,3) = ry;


% Stretch colormap according to minimum and maximum values
X = linspace(minval,maxval,length(cmap));

col = zeros(length(val),3);
for i = 1 : length(val)
    col(i,1) = interp1(X,cmap(:,1),val(i));
    col(i,2) = interp1(X,cmap(:,2),val(i));
    col(i,3) = interp1(X,cmap(:,3),val(i));
end

end
    
    






