%clear all

notload=false;

fdir='/Volumes/DISK_2020_5/ROCKY_BEACH/lloyd_lg_domain/';
m=656;
n=200;
l=30;


dx=0.015;
dy=0.015;
dz=0.0018;

x=[0:m-1]*dx;
y=[0:n-1]*dy;
z=-0.054+[0:l-1]*dz;

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


if(notload)

dep=load('depth_656x200.txt');

fnum=sprintf('%.5d',nfile(kt));

u=load([fdir 'u_' fnum]);
v=load([fdir 'v_' fnum]);
w=load([fdir 'w_' fnum]);
eta=load([fdir 'eta_' fnum]);
u3d1=reshape(u,[n,l,m]);
v3d1=reshape(v,[n,l,m]);
w3d1=reshape(w,[n,l,m]);

U3d=permute(u3d1,[1 3 2]);
V3d=permute(v3d1,[1 3 2]);
W3d=permute(w3d1,[1 3 2]);

u2d=squeeze(u3d1(:,l,:));
v2d=squeeze(v3d1(:,l,:));

end

sky=6;
skz=4;

nyz=[50 300 500];
for kkk=1:3
n1=nyz(kkk);
u_yz{kkk}=squeeze(U3d(1:sky:end,n1,1:skz:end));
v_yz{kkk}=squeeze(V3d(1:sky:end,n1,1:skz:end));
w_yz{kkk}=squeeze(W3d(1:sky:end,n1,1:skz:end));
x_yz{kkk}=squeeze(XXX(1:sky:end,n1,1:skz:end));
y_yz{kkk}=squeeze(YYY(1:sky:end,n1,1:skz:end));
z_yz{kkk}=squeeze(ZZZ(1:sky:end,n1,1:skz:end));
end


[vort vort1]=curl(X2d,Y2d,u2d,v2d);


%surf(X2d,Y2d,eta); shading interp;
   %alpha(0.1)
    colormap(bone)
    light('position',[-5,-5,0],'color','w');
    lighting phong
    material([0.8,0.8,0.5,10,0.5])
    hold on

surf(X2d,Y2d,-dep,'edgecolor','none','FaceColor',[0.9 1.0 1.0]);

%surfl(X2d,Y2d,-dep),shading interp;
   % alpha(0.7)
    colormap(bone)


sc=1.5;
scx=2.0;

for kkk=1:3
quiver3(x_yz{kkk},y_yz{kkk},z_yz{kkk},u_yz{kkk}*scx,v_yz{kkk}*sc,w_yz{kkk}*sc,0,'r','linewidth',1);
x1=x_yz{kkk}(1,1);
x2=x_yz{kkk}(end,1);
x3=x_yz{kkk}(end,end);
x4=x_yz{kkk}(1,end);
y1=y_yz{kkk}(1,1);
y2=y_yz{kkk}(end,1);
y3=y_yz{kkk}(end,end);
y4=y_yz{kkk}(1,end);
z1=z_yz{kkk}(1,1);
z2=z_yz{kkk}(end,1);
z3=z_yz{kkk}(end,end);
z4=z_yz{kkk}(1,end);

if(kkk==1)
plot3([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[z1 z2 z3 z4 z1],'k','LineWidth',1)
else
plot3([x2 x3 x4 x1],[y2 y3 y4 y1],[ z2 z3 z4 z1],'k','LineWidth',1)
end
end

% frame, n,m,k
x1=XXX(1,1,end);
x2=XXX(1,end,end);
x3=XXX(end,end,end);
x4=XXX(end,1,end);
y1=YYY(1,1,end);
y2=YYY(1,end,end);
y3=YYY(end,end,end);
y4=YYY(end,1,end);
zz=0.005;
z1=ZZZ(1,1,end)+zz;
z2=ZZZ(1,end,end)+zz;
z3=ZZZ(end,end,end)+zz;
z4=ZZZ(end,1,end)+zz;

%plot3([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[z1 z2 z3 z4 z1],'k-','LineWidth',1)
%    alpha(0.5)

vw=[-14 49];
view(vw);

grid
xlabel('x (m) ')
ylabel('y (m) ')
zlabel('z (m) ')

fname=['./plots/' 'flow3d_' fnum];

print('-djpeg100', fname)
print('-depsc2', [fname '.eps'] )












