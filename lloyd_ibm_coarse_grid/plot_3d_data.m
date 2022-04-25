clear all
fdir='/Users/fengyanshi/TMP/tmp1/';

m=164;
n=25;
l=30;

dx=0.06;
dy=0.06;
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
clf

for kt=1:length(nfile);

fnum=sprintf('%.5d',nfile(kt));

u=load([fdir 'u_' fnum]);
v=load([fdir 'v_' fnum]);
w=load([fdir 'w_' fnum]);
eta=load([fdir 'eta_' fnum]);
u3d1=reshape(u,[n,l,m]);
v3d1=reshape(v,[n,l,m]);
w3d1=reshape(w,[n,l,m]);

u3d=permute(u3d1,[1 3 2]);
v3d=permute(v3d1,[1 3 2]);
w3d=permute(w3d1,[1 3 2]);

u2d=squeeze(u3d1(:,l,:));
v2d=squeeze(v3d1(:,l,:));
[vort vort1]=curl(X2d,Y2d,u2d,v2d);

pcolor(X2d,Y2d,vort),shading flat
colorbar
caxis([-1 1])
title(['time= ' num2str(nfile(kt)) ' s'])
end

%print('-djpeg100', './plots/300s_500s_53cm.jpg')










