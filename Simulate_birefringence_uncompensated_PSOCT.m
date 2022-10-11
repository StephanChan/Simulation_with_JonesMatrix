%% Jones matrix calculation for PSOCT system
%x is the retardation, t is the fast axis orientation, 
syms x t;
syms y R sx st;
syms k zs zr;
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

for d=1:18
    d
%% define Jones matrix for optics
J_QWP_ref=subs(Jones,x,pi/2);
J_QWP_ref=subs(J_QWP_ref,t,pi/8);
J_samp_arm=subs(Jones,x,70/180*pi);
J_samp_arm=subs(J_samp_arm,t,10/180*pi);
J_QWP_samp=subs(Jones,x,pi/2);
J_QWP_samp=subs(J_QWP_samp,t,pi/4);

J_samp=subs(Jones,x,(25/180*pi));
J_samp=subs(J_samp,t,(d*10-10)/180*pi);
%% reference arm, zr=0
E_ref=simplify(J_QWP_ref*normal_glass*J_QWP_ref*[1;0]/2);

%% sample arm,zs=0
%signal_samp=simplify(J_QWP2*Jones*Jones*J_QWP2*[1;0]*sqrt(R)/2);
E_samp=simplify(J_QWP_samp*J_samp_arm*J_samp*J_samp*J_samp_arm*J_QWP_samp*[1;0]/2);

%% interference
channel1=simplify(2*E_ref(1)*conj(E_samp(1)),'Steps',50);
channel2=2*E_ref(2)*conj(E_samp(2));

%% retardance
ret(d)=vpa(atan(abs(channel1)/abs(channel2))/pi*180);
ori(d)=vpa(phase(channel1)-phase(channel2))/pi*180;
end
ori2=ori;
ori2(ori2<0)=ori2(ori2<0)+180;
ori2=ori2/2;
ori2(ori2<0)=ori2(ori2<0)+90;
ori2(ori2>45)=ori2(ori2>45)-90;
figure;plot((0:d-1)*10,ori2)
figure;plot((0:d-1)*10,ret)
%% sample arm with extra retarder, glass sample
%extra retarder is QWP, fast axis is 37 degree with respect to horizontal
% %plane
% J_QWP_extra=subs(Jones,x,pi/2);
% J_QWP_extra=subs(J_QWP_extra, t, 30/180*pi);
% E_samp=simplify(J_QWP_samp*normal_glass*J_QWP_samp*[0;1])

%% simulate cross polarization extinction ratio
% ratio=zeros(1,180);
% for i=1:180
%     J_QWP_samp=subs(Jones,x,pi/2);
%     J_QWP_samp=subs(J_QWP_samp,t,i/180*pi);
%     E_samp=simplify(J_QWP_samp*normal_glass*J_QWP_samp*[1;0]);
%     r=double(abs(E_samp(1))/abs(E_samp(2)));
%     ratio(i)=r;
% end
% figure;plot(abs(ratio));

%% animation of polarization
% t=0:0.005:1;
% x=sin(2*pi*t*5);
% y=0*sin(2*pi*t*2.31);
% % curve=animatedline('Marker','o');
% % set(gca);
% % grid on;
% % for i=1:length(t)
% %     addpoints(curve,x(i),y(i));
% %     drawnow
% % end
% p0=[0,0];
