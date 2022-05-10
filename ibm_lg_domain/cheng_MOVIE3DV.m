clear
clf
%%%%%%%%%%%%%%%%%%%% Use this code to make 3d movie %%%%%%%%%%%%%%%%%%%%%
ETA0=load('eta_03501');
DEP=-load('depth');
U=load('u_03501');
V=load('v_03501');
W=load('w_03501');
probe1=load('probe_0001');
probe2=load('probe_0002');
%%%%%%inputs%%%%%%
OUTPUT='VEC';%VEC(velocity field) or TUB(vortex tube)
D=0.054;%max wave depth
H=0.01;%max height
T=470;%initial time
TL=480;%ending time
DELT=0.1;%time interval
DELX=0.015;%spatial interval in X direction
DELY=0.015;%spatial interval in Y direction
NZ=30;%layer numbers in Z direction
LAYER=15;%layer numbers in Z direction for plotting
OBSX1=401;%the point 1 in x we observed the transect
OBSX2=367;%the point 1 in x we observed the transect
GAUGE1=51;%the point of gauge1 in y
GAUGE2=69;%the point of gauge1 in y
NX=length(ETA0(1,:));
NY=length(ETA0(:,1));
XL=NX*DELX;%spatial limitation in X direction
YL=NY*DELY;%spatial limitation in Y direction
%%%%%%%%% drawing variables %%%%%%%%%%%%%%
XI=DELX/2:DELX:XL-DELX/2;
YI=DELY/2:DELY:YL-DELY/2;
[XG,YG]=meshgrid(XI,YI);
XVEC1=repmat(XI(OBSX1)*ones(1,NY),LAYER,1);
XVEC2=repmat(XI(OBSX2)*ones(1,NY),LAYER,1);
YVEC=repmat(YI,LAYER,1);
%%%%%%%%% movie variables %%%%%%%%%%%%%%%%%
FDELT=DELT;%time interval for movie
DPF=FDELT/DELT;%using one frame per DPF steps
FPS=1/FDELT;%frame per second
% WriterObj=VideoWriter('SUBMERGED_3DVEL.avi');
WriterObj=VideoWriter('SUBMERGED_3DVCEL.avi');
WriterObj.FrameRate=FPS;
open(WriterObj);
%%%%%%%%%%%%%%%%%%start making%%%%%%%%%%%%%%
NUM=TL/DELT;
[ICTX,ICTY]=find(abs(ETA0)>10e-4);%find where is land
ICT=[ICTX,ICTY];
ICTN=length(ICT);
for n=T/DELT+1:NUM
  FileName=['eta_' num2str(n,'%05d')];
  ETA=load(FileName);
  
  if OUTPUT=='VEC'
  %%%%%%%%%%%%%3D Velocity vectors
  %%%%%%%X,Y,Z,U,V,W at point1%%%%%%%%%
  ZVEC1=zeros(LAYER,NY);
  for i=1:NY
    DELZ=(ETA(i,OBSX1)-DEP(i,OBSX1))/LAYER;
    ZVEC1(:,i)=linspace(DEP(i,OBSX1)+DELZ/2,ETA(i,OBSX1)-DELZ/2,LAYER);
  end
  FileName=['u_' num2str(n,'%05d')];
  U=load(FileName);
  u3d=reshape(U,[NY,NZ,NX]);
  UVEC1=squeeze(u3d(:,2:2:end,OBSX1));
  UVEC1=UVEC1';
  FileName=['v_' num2str(n,'%05d')];
  V=load(FileName);
  v3d=reshape(V,[NY,NZ,NX]);
  VVEC1=squeeze(v3d(:,2:2:end,OBSX1));
  VVEC1=VVEC1';
  FileName=['w_' num2str(n,'%05d')];
  W=load(FileName);
  w3d=reshape(W,[NY,NZ,NX]);
  WVEC1=squeeze(w3d(:,2:2:end,OBSX1));
  WVEC1=WVEC1';
  %%%%%%%X,Y,Z,U,V,W at point2%%%%%%%%%
  ZVEC2=zeros(LAYER,NY);
  for i=1:NY
    DELZ=(ETA(i,OBSX2)-DEP(i,OBSX2))/LAYER;
    ZVEC2(:,i)=linspace(DEP(i,OBSX2)+DELZ/2,ETA(i,OBSX2)-DELZ/2,LAYER);
  end
  UVEC2=squeeze(u3d(:,2:2:end,OBSX2));
  UVEC2=UVEC2';
  VVEC2=squeeze(v3d(:,2:2:end,OBSX2));
  VVEC2=VVEC2';
  WVEC2=squeeze(w3d(:,2:2:end,OBSX2));
  WVEC2=WVEC2';
  %%%%%%%set NaN for land%%%%%%%%%
  for k=1:ICTN
    if ETA(ICT(k,1),ICT(k,2))==ETA0(ICT(k,1),ICT(k,2))
      ETA(ICT(k,1),ICT(k,2))=NaN;%set NaN for land
    else
    end
  end  
  if rem(n,DPF)==0
    TIME=n*DELT;
    %%%%%%% draw picture %%%%%%%%%%%%%%%%%
    h = figure(1);
    pp=[9 9 1904 1041];set(gcf,'position',pp)
%     set(gcf,'outerposition',get(0,'screensize'));
    set(h,'Color','white')
%     %X-Z velocity
%     f2=subplot(2,2,2);
%     UWM=max(sqrt(UVEC1(:,NY/2).^2+WVEC1(:,NY/2).^2));
%     set(f2,'Units', 'normalized', 'Position', [0.8 0.55 0.15 0.4])%[left£¬bottom£¬width£¬height]
%     set(f2,'fontsize',20,'linewidth',2,'box','on')
%     qscale = 0.35 ; % scaling factor for all vectors
%     h2=quiver(XVEC1(:,NY/2),ZVEC1(:,NY/2),UVEC1(:,NY/2),WVEC1(:,NY/2),0,'r','linewidth',2);
%     hU = get(h2,'UData') ;
%     hV = get(h2,'VData') ;
%     set(h2,'UData',qscale*hU,'VData',qscale*hV)
%     xlabel('X','fontsize',15)
%     ylabel('Z','fontsize',15)
%     uwm=['Below Gauge 1'  '            max |U|: ' num2str(UWM,'%.4f') ' (m/s)'];
%     title(uwm)
%     xlim([XI(OBSX1)-0.01,XI(OBSX1)+0.04])
%     ylim([-0.054 0.01])
%     set(f2,'fontsize',12,'linewidth',2,'box','on')
    %%%stacked uw plane
    f2=subplot(2,2,2);
    UWM=max(sqrt(UVEC1(:,NY/2).^2+WVEC1(:,NY/2).^2));
    set(f2,'Units', 'normalized', 'Position', [0.8 0.55 0.15 0.4])%[left£¬bottom£¬width£¬height]
    set(f2,'fontsize',20,'linewidth',2,'box','on')
    X1=zeros(LAYER,1);
    Z1=[0.01:0.01:0.15]';
    h2=quiver(X1,Z1,UVEC1(:,NY/2),WVEC1(:,NY/2),0,'r','linewidth',2);
    xlabel('u(m/s)','fontsize',15)
    ylabel('w(m/s)','fontsize',15)
    uwm=['Below Gauge 1'  '            max |U|: ' num2str(UWM,'%.4f') ' (m/s)'];
    title(uwm)
    xlim([-0.01,0.13])
    ylim([0 0.16])
    set(f2,'fontsize',12,'linewidth',2,'box','on')
    %Y-Z velocity
%     f3=subplot(2,2,4);
%     VWM=max(sqrt(VVEC1(:,NY/2).^2+WVEC1(:,NY/2).^2));
%     set(f3,'Units', 'normalized', 'Position', [0.8 0.05 0.15 0.4])%[left£¬bottom£¬width£¬height]
%     set(f3,'fontsize',20,'linewidth',2,'box','on')
%     h3=quiver(YVEC(:,NY/2),ZVEC1(:,NY/2),VVEC1(:,NY/2),WVEC1(:,NY/2),0,'r','linewidth',2);
%     hU = get(h3,'UData') ;
%     hV = get(h3,'VData') ;
%     set(h3,'UData',qscale*hU,'VData',qscale*hV)
%     xlabel('Y','fontsize',15)
%     ylabel('Z','fontsize',15)
%     vwm=['                         max |U|: ' num2str(VWM,'%.4f') ' (m/s)'];
%     title(vwm)
%     xlim([YI(NY/2)-0.025,YI(NY/2)+0.025])
%     ylim([-0.054 0.01])
%     set(f3,'fontsize',12,'linewidth',2,'box','on')
    %%%stacked vw plane
    f3=subplot(2,2,4);
    VWM=max(sqrt(VVEC1(:,NY/2).^2+WVEC1(:,NY/2).^2));
    set(f3,'Units', 'normalized', 'Position', [0.8 0.05 0.15 0.4])%[left£¬bottom£¬width£¬height]
    set(f3,'fontsize',20,'linewidth',2,'box','on')
    Y2=zeros(LAYER,1);
    Z2=[0.01:0.01:0.15]';
    h3=quiver(Y2,Z2,VVEC1(:,NY/2),WVEC1(:,NY/2),0,'r','linewidth',2);
    xlabel('v(m/s)','fontsize',15)
    ylabel('w(m/s)','fontsize',15)
    vwm=['                         max |U|: ' num2str(VWM,'%.4f') ' (m/s)'];
    title(vwm)
    xlim([-0.08,0.08])
    ylim([0 0.16])
    set(f3,'fontsize',12,'linewidth',2,'box','on')
    %3D velocity
    f1=subplot(2,2,[1,3]);
    set(f1,'Units', 'normalized', 'Position', [0.05 0.1 0.7 0.8])%[left£¬bottom£¬width£¬height]
    set(f1,'fontsize',20,'linewidth',2,'box','on')
    set(f1,'DataAspectRatio',[1 1 0.1])
    axis([4 7 0 YL -D H]);
    view(-25,30)
    hold on
    surf(XG,YG,ETA); shading interp;
    alpha(0.3)
    colormap(bone)
    light('position',[-5,-5,5],'color','w');
    lighting phong
    material([0.8,0.8,0.5,10,0.5])
    hold on
    surf(XG,YG,DEP,'edgecolor','none','FaceColor',[0.8 1.0 0.8]);
    plot([0 10],[2 2],'linewidth',2);
    xlabel('$\begin{array}{c} X \end{array}$','interpreter','latex','fontsize',20)
    ylabel('$\begin{array}{c} Y \end{array}$','interpreter','latex','fontsize',20)
    zlabel('$\begin{array}{c} Z \end{array}$','interpreter','latex','fontsize',20)
    set(f1,'fontsize',20,'linewidth',2,'box','on','BoxStyle','full','Color',[0.8 1.0 0.8]);
    %point 1
    %%%%%%% draw 3d vectors %%%%%%%%%%%%%%
    quiver3(XVEC1(1:2:end,1:5:end),YVEC(1:2:end,1:5:end),ZVEC1(1:2:end,1:5:end),UVEC1(1:2:end,1:5:end),VVEC1(1:2:end,1:5:end),WVEC1(1:2:end,1:5:end),0.8,'r','linewidth',2);
    plot3([XI(OBSX1) XI(OBSX1)],[0 0],[-D H],'k','linewidth',2);
    plot3([XI(OBSX1) XI(OBSX1)],[YL YL],[-D H],'k','linewidth',2);
    plot3([XI(OBSX1) XI(OBSX1)],[0 YL],[-D -D],'k','linewidth',2);
    plot3([XI(OBSX1) XI(OBSX1)],[0 YL],[H H],'k','linewidth',2);
    %%%%%%% mark two gauges %%%%%%%%%%%%
    plot3(XVEC1(LAYER,GAUGE1),YVEC(LAYER,GAUGE1),ETA(GAUGE1,OBSX1),'o','Markerfacecolor','k','markeredgecolor','k')%gauge 1
    text(XVEC1(LAYER,GAUGE1),YVEC(LAYER,GAUGE1),ETA(GAUGE1,OBSX1)+DELZ,'G1','fontsize',15)
    plot3(XVEC1(LAYER,GAUGE1),YVEC(LAYER,GAUGE2),ETA(GAUGE2,OBSX1),'o','Markerfacecolor','k','markeredgecolor','k')%gauge 2
    text(XVEC1(LAYER,GAUGE1),YVEC(LAYER,GAUGE2),ETA(GAUGE2,OBSX1)+DELZ,'G2','fontsize',15)
    plot3(XVEC1(1,:),YI,ETA(:,OBSX1),'-b','linewidth',2)
    %point 2
    %%%%%%% draw 3d vectors %%%%%%%%%%%%%%
    quiver3(XVEC2(1:2:end,1:5:end),YVEC(1:2:end,1:5:end),ZVEC2(1:2:end,1:5:end),UVEC2(1:2:end,1:5:end),VVEC2(1:2:end,1:5:end),WVEC2(1:2:end,1:5:end),0.8,'r','linewidth',2);
    plot3([XI(OBSX2) XI(OBSX2)],[0 0],[-D H],'k','linewidth',2);
    plot3([XI(OBSX2) XI(OBSX2)],[YL YL],[-D H],'k','linewidth',2);
    plot3([XI(OBSX2) XI(OBSX2)],[0 YL],[-D -D],'k','linewidth',2);
    plot3([XI(OBSX2) XI(OBSX2)],[0 YL],[H H],'k','linewidth',2);
    %%%%%%% draw frame %%%%%%%%%%%%
    plot3(XVEC2(1,:),YI,ETA(:,OBSX2),'-b','linewidth',2)
    %%%%%%% draw title with time %%%%%%%%%%%%%%%%%
    tm=['TIME: ' num2str(TIME,'%.1f') ' (second)'];
    title(tm)
    %%%%%%% add time interval %%%%%%%%%%%%
    T=T+DELT;
    %%%%%%% save frame %%%%%%%%%%%%%%%%%%%
    FRAME(n/DPF) = getframe(gcf);
    writeVideo(WriterObj,FRAME(n/DPF))
    clf
  end
   
  elseif OUTPUT=='TUB'

  %%%%%%%%%%%%%3D vortex tube
  %obtain velocity
  FileName=['u_' num2str(n,'%05d')];
  U=load(FileName);
  u3d=reshape(U,[NY,NZ,NX]);
  U3D=permute(u3d,[1 3 2]);
  FileName=['v_' num2str(n,'%05d')];
  V=load(FileName);
  v3d=reshape(V,[NY,NZ,NX]);
  V3D=permute(v3d,[1 3 2]);
  FileName=['w_' num2str(n,'%05d')];
  W=load(FileName);
  w3d=reshape(W,[NY,NZ,NX]);
  W3D=permute(w3d,[1 3 2]);
  %location of velocty
  DELZ1=(ETA-DEP)/NZ;DELZ2(:,1,:)=DELZ1;DELZ=repmat(DELZ2,1,NZ);
  DEP1(:,1,:)=DEP;DEP2=repmat(DEP1,1,NZ);
  ZI=cumsum(DELZ,2)-DELZ/2+DEP2;
  DELZ=D/NZ;
  ZII=-D+DELZ/2:DELZ:-DELZ/2;
  [X3D,Y3D,Z3D2]=meshgrid(XI,YI,ZII);
  Z3D=permute(ZI,[1 3 2]);
  %calculate curl
  [curlx,curly,~,~] = curl(X3D,Y3D,Z3D,U3D,V3D,W3D);  
  if rem(n,DPF)==0
    TIME=n*DELT;
    %%%%%%% draw picture %%%%%%%%%%%%%%%%%
    h = figure(1);
    pp=[9 9 1904 1041];set(gcf,'position',pp)
%     set(gcf,'outerposition',get(0,'screensize'));
    set(h,'Color','white')
    slice(X3D,Y3D,Z3D2,curlx,[6 6.5 7],1.5,-D);shading interp
    colormap(b2r(-0.2,0.2));
%     slice(X3D,Y3D,Z3D2,curly,10,[0.75],-D);shading interp
%     colormap(b2r(-1.5,1.5));
    colorbar
    material([0.8,0.8,0.5,10,0.5])
    hold on
    surf(XG,YG,DEP,'edgecolor','none','FaceColor',[0.8 1.0 0.8]);
    p1 = patch(isosurface(X3D,Y3D,Z3D,curlx,0.2));
%     p1 = patch(isosurface(X3D,Y3D,Z3D,curly,1.5));
    alpha(p1,1)
    xlim([4,8])
    p1.FaceColor = 'red';
    p1.EdgeColor = 'none';
    p2 = patch(isosurface(X3D,Y3D,Z3D,curlx,-0.2));
%     p2 = patch(isosurface(X3D,Y3D,Z3D,curly,-1.5));
    alpha(p2,1)
    xlim([4,8])
    zlim([-D,0])
    p2.FaceColor = 'blue';
    p2.EdgeColor = 'none';
    daspect([1 1 0.1])
    view(3); 
    camlight 
    lighting gouraud
    xlabel('$x$(m)', 'interpreter', 'latex', 'fontsize',20)
    ylabel('$y$(m)', 'interpreter', 'latex', 'fontsize',20)
    zlabel('$z$(m)', 'interpreter', 'latex', 'fontsize',20)
    set(gca,'YTickLabel',{'-0.75' '-0.25' '0.25' '0.75'})
    set(gca,'fontname','Times','fontsize',20,'linewidth',2,'box','on')
    %%%%%%% draw title with time %%%%%%%%%%%%%%%%%
    tm=['TIME: ' num2str(TIME,'%.1f') ' (second)'];
    title(tm)
    %%%%%%% add time interval %%%%%%%%%%%%
    T=T+DELT;
    %%%%%%% save frame %%%%%%%%%%%%%%%%%%%
    FRAME(n/DPF) = getframe(gcf);
    writeVideo(WriterObj,FRAME(n/DPF))
    clf
  end
  
  end
end
close(WriterObj);