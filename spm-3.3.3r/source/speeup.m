function speedup()

close all;

%first test
test{1}.n = [1,4,12,24,48];
test{1}.t = [150,38.3,13.4,7.1,4.0];
test{1}.speedup = test{1}.t(1)./test{1}.t;

hold on; grid on; box on;
plot(test{1}.n,test{1}.speedup,'-ok');