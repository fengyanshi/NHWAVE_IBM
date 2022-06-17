clear all
m=656;
n=200;
k=30;
dep_flat=0.054;
%rat=1.1; % first test, no instability
rat=1.02;
dep(1:n,1:m)=dep_flat;

save -ASCII depth.txt dep


xstart=150;
ystart=20.0;
widthx=350.0;
widthy=100.0;

sizes=[60 70 80];
nums=[6 6 6];

rng(1);
for kk=1:length(sizes)
x0{kk}=floor(xstart+rand(nums(kk),1)*widthx);
x1{kk}=x0{kk}+sizes(kk);
y0{kk}=floor(ystart+rand(nums(kk),1)*widthy);
y1{kk}=y0{kk}+sizes(kk);
end

sum_zzz=zeros([n,m]);

figure(1)
wid=10;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf
hold on
for kk=1:length(sizes)
for kkk=1:nums(kk)
zzz=zeros([n,m]);
xx1=x0{kk}(kkk);
xx2=x1{kk}(kkk);
yy1=y0{kk}(kkk);
yy2=y1{kk}(kkk);
plot([xx1 xx2 xx2 xx1 xx1],[yy1 yy1 yy2 yy2 yy1]);
zz=peaks(sizes(kk));
zzz(yy1:yy1+sizes(kk)-1, xx1:xx1+sizes(kk)-1)=zz(:,:);
sum_zzz=sum_zzz+zzz;
end
end

sum_zzz(sum_zzz<0)=0.0;
% normalize 
mm=max(max(sum_zzz));
mn=min(min(sum_zzz));
sum_zzz=(sum_zzz-mn)/(mm-mn);

dep=dep-sum_zzz*dep_flat/rat;

print('-djpeg100', './plots/blocks.jpg')

figure(2)
wid=10;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf
colormap jet
subplot(211)
pcolor(sum_zzz),shading flat
title('normalized topography')
ylabel('points in y')
colorbar
subplot(212)
pcolor(-dep),shading flat
colorbar
title('bathymetry')
ylabel('points in y')
ylabel('points in x')

print('-djpeg100', './plots/depth_pcolor.jpg')

figure(3)
wid=10;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

mesh(-dep)
view([-18 60])
title('bathymetry')
ylabel('points in y')
xlabel('points in x')


print('-djpeg100', './plots/depth_mesh.jpg')

save -ASCII depth_656x200.txt dep


