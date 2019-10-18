% Analyzing multiple nonlinear time series with extended Granger causality
% Examples
addpath('/home/xyy/matcode/GC_clean/GCcal/');

% Example 1. Eq. (11)

n = 1e6;

% "In order to make the simulations realistic, some system noise and measurement noise are added to the time" series.

x = zeros(n, 1);
y = zeros(n, 1);
%err = 1e-6 * rand(n, 2);

err0 = 0.3;
x(1) = err0 * rand(1);
x(2) = err0 * rand(1);
y(1) = err0 * rand(1);
y(2) = err0 * rand(1);

c = 0.5;

for k = 3:n
  x(k) = 3.4 * x(k-1) * (1 - x(k-1)^2) * exp(-x(k-1)^2) + 0.8 * x(k-2);
  y(k) = 3.4 * y(k-1) * (1 - y(k-1)^2) * exp(-y(k-1)^2) + 0.5 * y(k-2) + c * x(k-2);
end

od = 20;

GC = nGrangerTfast([x'; y'], od);
disp('GC ex1.1')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex1.1 nonzero prob')
disp(nzp)

% Pic output
figure(1001);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex1.1', p_value, gc0);

figure(93);
rg = 2:1e5;
plot(x(rg), x(rg-1), '.');
xlabel('x(t)');
ylabel('x(t-1)');

figure(94);
rg = 2:1e5;
plot(y(rg), y(rg-1), '.');
xlabel('y(t)');
ylabel('y(t-1)');

% Example 1. Eq. (12)
for k = 3:n
  x(k) = 3.4 * x(k-1) * (1 - x(k-1)^2) * exp(-x(k-1)^2) + 0.8 * x(k-2);
  y(k) = 3.4 * y(k-1) * (1 - y(k-1)^2) * exp(-y(k-1)^2) + 0.5 * y(k-2) + c * x(k-2)^2;
end

od = 20;

GC = nGrangerTfast([x'; y'], od);
disp('GC ex1.2')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex1.2 nonzero prob')
disp(nzp)

% Pic output
figure(1002);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex1.2', p_value, gc0);

figure(95);
rg = 2:1e5;
plot(y(rg), y(rg-1), '.');
xlabel('y(t)');
ylabel('y(t-1)');



