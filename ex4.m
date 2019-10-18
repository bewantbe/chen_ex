% Analyzing multiple nonlinear time series with extended Granger causality
% Examples
addpath('/home/xyy/matcode/GC_clean/GCcal/');
pic_output = @(st) print('-depsc2', sprintf('%s.eps', st));

len = 1e5;
x = zeros(1, len);
y = zeros(1, len);
z = zeros(1, len);

err = randn(3, len);

for k = 2 : len
  x(k) = 0.2 * x(k-1) + err(1, k);
  y(k) = 0.5 * y(k-1) + 0.5 * x(k-1) + err(2, k);
  z(k) = 0.4 * z(k-1) + 0.3 * y(k-1) + c * x(k-1) + err(3, k);
end

od = 20;
GC = nGrangerTfast([x; y; z], od);
disp('GC ex4.1')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex4.1 nonzero prob')
disp(nzp)

% Pic output
figure(1005);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex4.1', p_value, gc0);

% ex4.2, Eq.(16)

x = zeros(1, len);
y = zeros(1, len);
z = zeros(1, len);

x(1) = 0.1 * randn(1);
y(1) = 0.1 * randn(1);
z(1) = 0.1 * randn(1);

for k = 2 : len
  x(k) = 3.4 * x(k-1) * (1 - x(k-1)^2) * exp(-x(k-1)^2);
  y(k) = 3.4 * y(k-1) * (1 - y(k-1)^2) * exp(-y(k-1)^2) + 0.5 * x(k-1);
  z(k) = 3.4 * z(k-1) * (1 - z(k-1)^2) * exp(-z(k-1)^2) + 0.3 * y(k-1) + c * x(k-1);
end

od = 20;
GC = nGrangerTfast([x; y; z], od);
disp('GC ex4.2')
disp(GC)

p = size(GC,1);
len = length(x);
nzp = gc_prob_nonzero(GC, od, len-(p+1)*od);
disp('GC ex4.2 nonzero prob')
disp(nzp)

% Pic output
figure(1006);
p_value = 0.001;
gc0 = chi2inv(1-p_value, od)/len;
SimpleGCPic(GC, 'ex4.2', p_value, gc0);

