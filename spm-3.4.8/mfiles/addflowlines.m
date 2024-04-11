function addflowlines()

load flines.mat
hold on;
for i=1:length(flines),
    plot3(flines{i}.x,flines{i}.y,flines{i}.z);
end