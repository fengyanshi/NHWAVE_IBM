clear all

Mloc=164;
Nloc=25;
Kloc=30;

Nghost=2;
Ibeg=1+Nghost;
Iend=Mloc-Nghost;
Jbeg=1+Nghost;
Jend=Nloc-Nghost;
Kbeg=1+Nghost;
Kend=Kloc-Nghost;

dx=0.06;
dy=0.06;
dz=0.54/(Kloc-1);

% get mask and rx and ry ------------------

   % default is fluid
    set_flag = zeros([Mloc,Nloc,Kloc]);
    D=zeros([Mloc,Nloc])+0.054;
    sigc=[0:Kloc-1]*1./(Kloc-1);
    xc=[0:Mloc-1]*dx+dx/2.0;
    yc=[0:Nloc-1]*dy+dy/2.0;


    % specify obstacle 
    for k = 1:Kloc
    for j = 1:Nloc
    for i = 1:Mloc
      zc = D(i,j)*sigc(k);


      % obstacle  x<=4.5 & y<=13.4
	  dist = sqrt((xc(i)-5.0)^2+(yc(j)-0.75)^2);
      if(dist<0.375) 
	    zobst = min(0.049,0.049*(1.0-(dist-0.025)/(0.375-0.025)));
		if(zc<zobst) 
                   set_flag(i,j,k) = 1;
                 end
      end

    end
    end
    end

    % distance from obstacle/fluid interface to neighboring fluid cell 
    set_dist_x = zeros([Mloc,Nloc,Kloc]) -1.0;    % +1.e+20;
    set_dist_y = zeros([Mloc,Nloc,Kloc]) -1.0;    % +1.e+20;
    set_dist_z = zeros([Mloc,Nloc,Kloc]) -1.0;    % +1.e+20;

    % find the distance
    for k = 1:Kloc
    for j = Jbeg:Jend
    for i = Ibeg:Iend

      % obstacle #1 (x-5.0)^2+(y-0.75)^2=0.1815^2
      % loop over all fluid cells
      if(set_flag(i,j,k)==0)
	      zc = D(i,j)*sigc(k);
	      dist = min(0.025,(0.049-zc)/0.049*(0.375-0.025)+0.025);
          % obstacle on the right
          if(set_flag(i+1,j,k)==1)
            xlin = 5.0-sqrt(dist^2-(yc(j)-0.75)^2);
            set_dist_x(i,j,k) = abs(xc(i)-xlin);
          end

          % obstacle on the left
          if(set_flag(i-1,j,k)==1)
            xlin = 5.0+sqrt(dist^2-(yc(j)-0.75)^2);
            set_dist_x(i,j,k) = abs(xc(i)-xlin);
          end

          % obstacle on the front
          if(set_flag(i,j-1,k)==1)
            ylin = 0.75+sqrt(dist^2-(xc(i)-5.0)^2);
            set_dist_y(i,j,k) = abs(yc(j)-ylin);
          end
		  
          % obstacle on the back                                                                                               
          if(set_flag(i,j+1,k)==1)
            ylin = 0.75-sqrt(dist^2-(xc(i)-5.0)^2);
            set_dist_y(i,j,k) = abs(yc(j)-ylin);
          end
		  
		  % obstacle under the fluid
        if(k>1)
          if(set_flag(i,j,k-1)==1)  
            dist = sqrt((xc(i)-5.0)^2+(yc(j)-0.75)^2);
            zlin = min(0.049,0.049*(1.0-(dist-0.025)/(0.375-0.025)));
            set_dist_z(i,j,k) = abs(zc-zlin);
          end
        end % end k>1
      end

    end
    end
    end

%    set_dist_z(set_dist_z>=0)=dz-set_dist_z(set_dist_z>=0);
%    set_dist_x(set_dist_x>=0)=dx-set_dist_x(set_dist_x>=0);
%    set_dist_y(set_dist_y>=0)=dy-set_dist_y(set_dist_y>=0);

mask1=permute(set_flag,[2 3 1]);
mask2=reshape(mask1,[],Mloc);

save -ASCII obs_mask.txt mask2

mask1=permute(set_dist_x,[2 3 1]);
mask2=reshape(mask1,[],Mloc);

save -ASCII obs_distx.txt mask2

mask1=permute(set_dist_y,[2 3 1]);
mask2=reshape(mask1,[],Mloc);

save -ASCII obs_disty.txt mask2

mask1=permute(set_dist_z,[2 3 1]);
mask2=reshape(mask1,[],Mloc);

save -ASCII obs_distz.txt mask2

% figure
figure(1)
ox=load('obs_distx.txt');
oy=load('obs_disty.txt');
oz=load('obs_distz.txt');
o_all=ox+oy+oz;
pcolor(o_all),shading flat
colorbar
title('Rx + Ry +Rz (positive for distances, negative for void)');
xlabel('M (grid points)')
ylabel('N x L')
print -djpeg100 check_Rxyz.jpg

figure(2)
om=load('obs_mask.txt');
pcolor(om),shading flat
colorbar
title('masks (1 for obstacle)');
xlabel('M (grid points)')
ylabel('N x L')
print -djpeg100 check_mask.jpg








