clear all
fdir='/Volumes/DISK_2020_5/ROCKY_BEACH/lloyd_lg_domain/';
m=656;
n=200;
l=30;

dep=load('depth_656x200.txt');

dx=0.015;
dy=0.015;
dz=0.0018;

x=[0:m-1]*dx;
y=[0:n-1]*dy;
z=[0:l-1]*dz;

[X2d Y2d]=meshgrid(x,y);

[X,Y,Z]=meshgrid(x,y,z);

nst=input('start= ');

nfile=[nst];
TIME={'300s','400s','500s'};

%nst=input('start=(155-220)');
%nen=input('end=(190)');

%f=figure(1);
wid=10;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf
colormap jet

for kt=1:length(nfile);

fnum=sprintf('%.5d',nfile(kt));

u=load([fdir 'u_' fnum]);
v=load([fdir 'v_' fnum]);
u3d1=reshape(u,[n,l,m]);
v3d1=reshape(v,[n,l,m]);

u3d=permute(u3d1,[1 3 2]);
v3d=permute(v3d1,[1 3 2]);

u2d=squeeze(u3d1(:,l,:));
v2d=squeeze(v3d1(:,l,:));
[vort vort1]=curl(X2d,Y2d,u2d,v2d);

pcolor(X2d,Y2d,vort),shading flat
colorbar
caxis([-1 1])
title(['time= ' num2str(nfile(kt)) ' s'])
hold on
vb=[-0.05 0.04 0.03 0.02 0.01];
contour(X2d,Y2d,-dep,'Color','k','LineWidth',1);
xlabel('x(m)')
ylabel('y(m)')

end

fname=['./plots/' 'vort_' fnum];

print('-djpeg100', fname)











