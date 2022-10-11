% %% Jones matrix calculation for PSOCT system
% %x is the retardation, t is the fast axis orientation, 
% syms x t;
% syms y R sx st;
% syms k zs zr;
% syms theta alpha;
% assume(theta,'real');
% assume(alpha,'real');
% assume(x,'real');
% assume(sx,'real');
% assume(t,'real');
% assume(st,'real');
% assume(zs,'real');
% assume(zr,'real');
% assume(k,'real');
% Jones=exp(-j*x/2)*[cos(t)^2+exp(j*x)*sin(t)^2 (1-exp(j*x))*cos(t)*sin(t);...
%     (1-exp(j*x))*cos(t)*sin(t) sin(t)^2+exp(j*x)*cos(t)^2];
% normal_glass=[-1 0; 0 -1];
% 
% LP=[cos(theta)^2,cos(theta)*sin(theta);cos(theta)*sin(theta),sin(theta)^2];
% %% define Jones matrix for optics
% 
% J_samp=subs(Jones,x,33/180*pi);
% J_samp=subs(J_samp,t,60/180*pi);
% J_QWP=subs(Jones,x,pi/2);
% J_QWP=subs(J_QWP,t,0/180*pi);
% LP=subs(LP,theta,45/180*pi);
% E_in_Jones=[cos(alpha);sin(alpha)];
% %% polarization propagation
% 
% Intensity=zeros(1,90);
% for ii=1:90
%     % rotate first linear polarizer
%     E_in=subs(E_in_Jones,alpha,ii/180*pi);
%     % polarization propagates
%     E_out=simplify(LP*J_QWP*J_samp*E_in);
%     % measure intensity
%     Intensity(ii)=abs(E_out(1))^2+abs(E_out(2))^2;
% end
% 
% %% fitting intensity to get retardance and phase
% func=@(p,zdata)p(1)*sin(p(2).*zdata./360*2*pi+p(3))+p(4);
% lb=[0,0,0,0];
% ub=[1,3,2*pi,1];
% opts = optimset('Display','off','TolFun',1e-10);
% 
% param=lsqcurvefit(func,[0.2,1,pi/2, 0.2],1:90,Intensity,lb,ub,opts);
% y=func(param,1:90);
% figure;
% scatter(1:90,Intensity);hold on;plot(1:90,y)
% xlabel('first linear polarizer orientation (degree)')
% ylabel('detected intensity')
% display(['retardance is: ',num2str(asin(param(1)/param(4))/pi*180), ' degree  orientation is: ', num2str(param(3)/pi*180/2),' degree'])
% 
% %% simulate cross polarization extinction ratio
% % ratio=zeros(1,180);
% % for i=1:180
% %     J_QWP_samp=subs(Jones,x,pi/2);
% %     J_QWP_samp=subs(J_QWP_samp,t,i/180*pi);
% %     E_samp=simplify(J_QWP_samp*normal_glass*J_QWP_samp*[1;0]);
% %     r=double(abs(E_samp(1))/abs(E_samp(2)));
% %     ratio(i)=r;
% % end
% % figure;plot(abs(ratio));
% 
% %% animation of polarization
% % t=0:0.005:1;
% % x=sin(2*pi*t*5);
% % y=0*sin(2*pi*t*2.31);
% % % curve=animatedline('Marker','o');
% % % set(gca);
% % % grid on;
% % % for i=1:length(t)
% % %     addpoints(curve,x(i),y(i));
% % %     drawnow
% % % end
% % p0=[0,0];

%% Simulate the linear polarizer experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Jones matrix calculation for PSOCT system
%x is the retardation, t is the fast axis orientation, 
syms x t;
syms y R sx st;
syms k zs zr;
syms theta alpha;
assume(theta,'real');
assume(alpha,'real');
assume(x,'real');
assume(sx,'real');
assume(t,'real');
assume(st,'real');
assume(zs,'real');
assume(zr,'real');
assume(k,'real');
Jones=exp(-j*x/2)*[cos(t)^2+exp(j*x)*sin(t)^2 (1-exp(j*x))*cos(t)*sin(t);...
    (1-exp(j*x))*cos(t)*sin(t) sin(t)^2+exp(j*x)*cos(t)^2];
normal_glass=[-1 0; 0 -1];

LP_Jones=[cos(theta)^2,cos(theta)*sin(theta);cos(theta)*sin(theta),sin(theta)^2];
%% define Jones matrix for optics
E_in_Jones=[cos(alpha);sin(alpha)];
E_in=subs(E_in_Jones,alpha,60/180*pi);

J_QWP=subs(Jones,x,pi/2);
J_QWP=subs(J_QWP,t,0/180*pi);
LP=subs(LP_Jones,theta,45/180*pi);

%% polarization propagation

Intensity=zeros(1,90);
for ii=1:90
    % rotate first linear polarizer
    J_samp=subs(LP_Jones,theta,ii/180*pi);
    % polarization propagates
    E_out=simplify(LP*J_QWP*J_samp*E_in);
    % measure intensity
    Intensity(ii)=abs(E_out(1))^2+abs(E_out(2))^2;
end

%% fitting intensity to get retardance and phase
func=@(p,zdata)p(1)*sin(p(2).*zdata./360*2*pi+p(3))+p(4);
lb=[0,0,0,0];
ub=[1,3,2*pi,1];
opts = optimset('Display','off','TolFun',1e-10);

param=lsqcurvefit(func,[0.2,1,pi/2, 0.2],1:90,Intensity,lb,ub,opts)
y=func(param,1:90);
figure;
scatter(1:90,Intensity);hold on;plot(1:90,y)
xlabel('sample linear polarizer orientation (degree)')
ylabel('detected intensity')
% display(['retardance is: ',num2str(asin(param(1)/param(4))/pi*180), ' degree  orientation is: ', num2str(param(3)/pi*180/2),' degree'])
