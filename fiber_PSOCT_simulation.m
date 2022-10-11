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

X=1:5:180;
V=zeros(length(X),2,2);
D=zeros(length(X),2);
jj=0;
for ii=1:5:180
    jj=jj+1
    %% define Jones matrix for optics
    J_in=subs(Jones,x,pi/9);
    J_in=subs(J_in,t,pi/8);
    J_samp=subs(Jones,x,ii/180*pi);
    J_samp=subs(J_samp,t,pi/8);
    J_surf=subs(Jones,x,pi/10);
    J_surf=subs(J_surf,t,pi/3);
    J_out=subs(Jones,x,pi/5);
    J_out=subs(J_out,t,pi/4);

    %% reference arm
    % E_ref=simplify(J_in*exp(-j*2*k*zr-j*pi)*normal_glass*J_in*[1;0]/2);

    %% sample arm
    %signal_samp=simplify(J_QWP2*Jones*Jones*J_QWP2*[1;0]*sqrt(R)/2);
    E_samp=J_out*J_samp*J_in*[1 0;0 1];
    E_surf=J_out*J_in*[1 0;0 1];
    %% remove J_in
    E=E_samp*E_surf';
    D(jj,:)=eig(E);

    %% retardance
    if imag(D(jj,1))>0
       ret(jj)=phase(D(jj,1)/D(jj,2))/2/pi*180;
    else
        ret(jj)=phase(D(jj,2)/D(jj,1))/2/pi*180;
    end
end
figure;plot(X,ret)
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
