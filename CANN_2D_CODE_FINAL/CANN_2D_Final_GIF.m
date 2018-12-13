%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-d Continuous Attractor Neural Network with Hebbian learning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
format long g
mm = 30; nn = 30;
%---------------------------------------------------------------%
Z = peaks;
surf(Z)
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,10) = 0;
k = 1;
%-------------------------------------------------------------%
%-----------------------------Weight matrix-------------------------------%
w_size = mm;
sig_w = 1.46;
w1=cell(w_size);
w=zeros(mm*nn);
count = 1;
C = 0.1;
Aw = 10;
dx = 2*pi/mm;
for loc_x=1:w_size;
    for loc_y=1:w_size;
        w1{loc_x,loc_y} =(Aw * weight_distance(loc_x,loc_y,w_size,w_size, dx,sig_w))-C;
        w(count,:) = reshape(w1{loc_x,loc_y}',1,mm*nn);
        count = count + 1;
    end
end
%-------------------------------------------------------------------------%

%---------------------------------External Input--------------------------%
sig_ext = 1;
I_ext=zeros(mm,nn);
A = 0.36;
I_ext = A * weight_distance(15,15,w_size,w_size, dx,sig_ext);
I = I_ext;
% for i = -3:3; 
% for j = -3:3;  
%      I_ext(mm/2+i,nn/2+j) = 10;
% end
% end
I_ext = reshape(I_ext',1,mm*nn);
%-------------------------------------------------------------------------%

%-----------------------u and r initialization----------------------------%
u0 = zeros(1,mm*nn);
u = u0;
r = zeros(1,mm*nn);
%-------------------------------------------------------------------------%

%----------------------------Divisive Norm of u---------------------------%
for i = 1:mm*nn
r = (u.^2) ./(1 + ( 0.5 * (sum(u.^2) - (u(i).^2)) ) );
end
% r = u - min(0, min(u(:)));
%-------------------------------------------------------------------------%

%--------------------------------updating u-------------------------------% 
tau_inv=1./10; 
dt = 1;
%---------------------------Iteration starts------------------------------%
for time=1:10000

%-----------------------External input is from 0 to 10--------------------%
if (time > 10)
   I_ext = zeros(1,mm*nn);
end

for i = 1:mm*nn
    u(i) =  tau_inv * dt * ((u(i)) + sum(w(i,:).*r) * dx + I_ext(i)) ;
end

%-------------------Divisive Norm of u afer updating it-------------------%
for m = 1:mm*nn
r = (u.^2) ./(1 + ( 0.5 * (sum(u.^2) - (u(m).^2)) ) );
end
% r = u - min(0, min(u(:)));

%-------------Plot the (activity) r after updating in surf----------------%
figure(1);
direction = [0 0 1];
h = surf(reshape(r',mm,nn));
coaa = colorbar('AxisLocation','out');
% get(aa,'Ytick');
% tt = [0:0.001:0.000002];
% set(aa,'Ytick',tt);
aa.Label.String = 'Activity Level';
caxis([0,0.006]);
view(0,90);
rotate(h,direction,180)
title(['Time=',num2str(time)]);
xlim([1,mm]);
ylim([1,nn]);
zlim([0,0.1]);
%--------------------------------------------------------------%
f = getframe(gcf);
im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
k = k+1;
%-------------------------------------------------------------%
pause(0);
end
%-------------------------------------------------------------------------%
imwrite(im,map,'CANN_2D_A_0.36.gif','DelayTime',0,'LoopCount',inf) %g443800