clear all

dep=load('depth_656x200.txt');
dep_flat=0.054;
dep=dep_flat-dep;

[n m]=size(dep);
k=30;
dx=0.015;
dy=0.015;
dz=dep_flat/(k-1);

z=[0:k-1]*dz;

set_flag=zeros([m n k]);
for kk=1:k
for i=1:m
for j=1:n
if z(kk)<dep(j,i)
set_flag(i,j,kk)=1;
end
end
end
end

mask1=permute(set_flag,[2 3 1]);
mask2=reshape(mask1,[],m);

save -ASCII obs_mask.txt mask2

% figure

figure(2)
om=load('obs_mask.txt');
pcolor(om),shading flat
colorbar
title('masks (1 for obstacle)');
xlabel('M (grid points)')
ylabel('N x L')
print -djpeg100 check_mask.jpg








