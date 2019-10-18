% Analyzing multiple nonlinear time series with extended Granger causality
% Examples
addpath('/home/xyy/matcode/GC_clean/GCcal/');
pic_output = @(st) print('-depsc2', sprintf('%s.eps', st));

c = 0.5;

odefunx = @(x1,y1,z1,x2,y2,z2) [
-(y1 + z1);
x1 + 0.2*y1;
0.2 + z1 * (x1 - 4.7);
-(y2 + z2) + c * x1;
x2 + 0.2 * y2;
0.2 + z2 * (x2 - 4.7) ];

odefun = @(y) odefunx(y(1),y(2),y(3),y(4),y(5),y(6));

y0 = randn(1, 6);

len = 10000;
yy = zeros(length(y0), len);

dt = 0.02;
yy(:, 1) = y0;
for j = 2:len
  yy(:, j) = yy(:, j-1) + dt * odefun(yy(:, j-1)) + sqrt(dt)*1e-1*rand(6,1);
end

x = yy(1, :);
y = yy(2, :);
z = yy(3, :);

od = 20;
GC = nGrangerTfast([x; y], od);
disp('GC ex3')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex3 nonzero prob')
disp(nzp)

% Pic output
figure(1004);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex3', p_value, gc0);

figure(131);
rg = 1:length(x);
plot(rg, [x(rg)', y(rg)']);

