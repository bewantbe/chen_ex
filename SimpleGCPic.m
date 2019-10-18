% simple GC pic output

function SimpleGCPic(GC, fname, p_value, gc0)

pic_output = @(st) print('-depsc2', sprintf('%s.eps', st));

p = length(GC);

if p == 2
  gcs = [GC(2,1), GC(1,2)]
  lbs = {'x->y', 'y->x'};
elseif p == 3
  gcs = GC(eye(3)==0);
  lbs = {'x->y', 'x->z', 'y->x', 'y->z', 'z->x', 'z->y'};
end

ngc = length(gcs);

bar(1:ngc, gcs);
hold on
ax = plot([0,ngc+1], ones(1,2) * gc0);
legend(ax, sprintf('p-val = %.1e', p_value));
set(gca, 'xtick', 1:ngc, 'xticklabel', lbs);
if max(gcs) > 10*gc0
  set(gca, 'yscale', 'log');
end
ylabel('GC');
hold off
pic_output(fname);

end % function
