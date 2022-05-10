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

[XXX,YYY,ZZZ]=meshgrid(x,y,z);

%nst=input('start= ');
nst=500;

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

kt=1;

fnum=sprintf('%.5d',nfile(kt));

u=load([fdir 'u_' fnum]);
v=load([fdir 'v_' fnum]);
w=load([fdir 'w_' fnum]);
u3d1=reshape(u,[n,l,m]);
v3d1=reshape(v,[n,l,m]);
w3d1=reshape(w,[n,l,m]);

U3d=permute(u3d1,[1 3 2]);
V3d=permute(v3d1,[1 3 2]);
W3d=permute(w3d1,[1 3 2]);

u2d=squeeze(u3d1(:,l,:));
v2d=squeeze(v3d1(:,l,:));
[vort vort1]=curl(X2d,Y2d,u2d,v2d);

%[dudx,dudy,dudz]=gradient(U3d,dx,dy,dz);
%[dvdx,dvdy,dvdz]=gradient(V3d,dx,dy,dz);
%[dwdx,dwdy,dwdz]=gradient(W3d,dx,dy,dz);

%q=-0.5*(dudx.^2+dvdy.^2+dwdz.^2)-dudy.*dvdx-dudz.*dwdx-dvdz.*dwdy;

vw=[36 90];

[faces,verts,colors]=isosurface(XXX,YYY,ZZZ,q,0.001,ZZZ);
patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'none')
view(vw)
camlight; lighting gouraud
grid
xlabel('x (m) ')
ylabel('y (m) ')


fname=['./plots/' 'q_' fnum];

%print('-djpeg100', fname)











