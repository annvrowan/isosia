function status()

%  DEMstatus() writes status of current model run
%
%    DLH - 04.05.05
%

st = load('./status.dat');

disp([' ']);
disp(['Status: time = ',num2str(st(1)),', ',num2str(st(2)),'% done, filenr = ',num2str(st(3))]); 
%eval(['type ./output/logfile.txt']);
disp([' ']);