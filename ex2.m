% Analyzing multiple nonlinear time series with extended Granger causality
% Examples
addpath('/home/xyy/matcode/GC_clean/GCcal/');
pic_output = @(st) print('-depsc2', sprintf('%s.eps', st));

c = 0.5;

odefun = @(t, y) [
  -0.25*y(1) + y(2) - y(2)^3;
  y(1) - y(2) - y(1)*y(2);
  -0.25*y(3) + y(4) - y(4)^3 + c * y(1);
  y(3) - y(4) - y(3)*y(4)];
y0 = [1.6130    0.6268    4.3849    0.8178];
%[tt, yy] = ode45(odefun, linspace(0, 5*1000, 100*1000), y0);

len = 1e5;
yy = zeros(4, len);
yy(:,1) = y0;
dt = 0.05;
for j=2:len
  yy(:,j) = yy(:,j-1) + dt * odefun(j, yy(:,j-1)) + sqrt(dt)*0.01*rand(4,1);
end

x = yy(1, :);
y = yy(3, :);

GC = nGrangerTfast([x; y], 20);
disp('GC ex2')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex3 nonzero prob')
disp(nzp)

% Pic output
figure(1003);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex2', p_value, gc0);

%figure(101);
%rg = 1:length(x);
%plot(rg, [x(rg)', y(rg)']);

%figure(102);
%rg = 1:length(x);
%plot(x(rg), y(rg), '.');

