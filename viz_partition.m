function viz_partition(N)

col = {'*', 'r*', 'k*', 'm*', 'g*', 'y*'};

figure; hold on;
for i=0:N-1
  file = ['/tmp/rafb/testoutput/TestDistributedQuadMeshPartitioning/res_',num2str(N),'_',num2str(i),'.txt'];
  d = load(file);
  plot(d(:,1),d(:,2),col{i+1});
end;

