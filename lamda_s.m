clear;
clc;
epsilon=0.6;
UpperDelta=1;
barK=0;
phi=max(epsilon,UpperDelta);
gamma=0.1;
N=5;
A=[0 0 0 0 1;
   1 0 1 1 0;
   0 1 0 0 0;
   0 0 1 0 0;
   1 0 1 0 0];

C=10.0843;

s=1;

cvx_begin gp
    variable lamda
    variables bI(5) 
    variables bE(5) 
    variables deltah(5)
    variables tao(5) 
    variables u(5)
    variable v(5)
    obj=0;

    for ii=1:N
         obj=obj+hI(bI(ii))+hE(bE(ii))+ghat(phi,deltah(ii))+fhat(gamma,tao(ii));
    end


    minimize(lamda-phi)
    subject to
        obj<=C*s;
        epsilon*u(1)+deltah(1)*v(1)<=lamda*v(1);
        epsilon*u(2)+deltah(2)*v(2)<=lamda*v(2);
        epsilon*u(3)+deltah(3)*v(3)<=lamda*v(3);
        epsilon*u(4)+deltah(4)*v(4)<=lamda*v(4);
        epsilon*u(5)+deltah(5)*v(5)<=lamda*v(5);
        tao(1)*(bE(1)*u(5)+bI(1)*v(5))+(phi-epsilon)*u(1)<=lamda*u(1);
        tao(2)*(bE(2)*u(1)+bI(2)*v(1))+tao(2)*(bE(2)*u(3)+bI(2)*v(3))+tao(2)*(bE(2)*u(4)+bI(2)*v(4))+(phi-epsilon)*u(2)<=lamda*u(2);
        tao(3)*(bE(3)*u(2)+bI(3)*v(2))+(phi-epsilon)*u(3)<=lamda*u(3);
        tao(4)*(bE(4)*u(4)+bI(4)*v(4))+(phi-epsilon)*u(4)<=lamda*u(4);
        tao(5)*(bE(5)*u(1)+bI(5)*v(1))+tao(5)*(bE(5)*u(3)+bI(5)*v(3))+(phi-epsilon)*u(5)<=lamda*u(5);
        
        gamma/(1+gamma)<=tao(1)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(2)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(3)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(4)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(5)<=gamma/(0.1+gamma);
        phi-1<=deltah(1)<=phi-0.1;
        phi-1<=deltah(2)<=phi-0.1;
        phi-1<=deltah(3)<=phi-0.1;
        phi-1<=deltah(4)<=phi-0.1;
        phi-1<=deltah(5)<=phi-0.1;
        0.06<=bI(1)<=0.6;
        0.06<=bI(2)<=0.6;
        0.06<=bI(3)<=0.6;
        0.06<=bI(4)<=0.6;
        0.06<=bI(5)<=0.6;
        0.1<=bE(1)<=0.7;
        0.1<=bE(2)<=0.7;
        0.1<=bE(3)<=0.7;
        0.1<=bE(4)<=0.7;
        0.1<=bE(5)<=0.7;        
cvx_end

delta=phi*ones(N,1)-deltah;
theta=gamma*(1-tao)./tao;



BI=diag(bI);
BE=diag(bE);
E=epsilon*eye(N);
D=diag(delta);
T=diag(gamma*ones(N,1)./(theta+gamma*ones(N,1)));
Q=[T*BE*A-E T*BI*A;
    E -D];

max(real(eig(Q)))


syms v1 v2 v3 v4 v5 e1 e2 e3 e4 e5 i1 i2 i3 i4 i5 s1 s2 s3 s4 s5
eq1=gamma*v1-theta(1)*s1-s1*(bE(1)*e5+bI(1)*i5)==0;
eq2=s1*(bE(1)*e5+bI(1)*i5)-epsilon*e1==0;
eq3=epsilon*e1-delta(1)*i1==0;
eq4=delta(1)*i1+theta(1)*s1-gamma*v1==0;
eq5=i1+v1+s1+e1==1;

eq6=gamma*v2-theta(2)*s2-s2*(bE(2)*e1+bI(2)*i1+bE(2)*e3+bI(2)*i3+bE(2)*e4+bI(2)*i4)==0;
eq7=s2*(bE(2)*e1+bI(2)*i1+bE(2)*e3+bI(2)*i3+bE(2)*e4+bI(2)*i4)-epsilon*e2==0;
eq8=epsilon*e2-delta(2)*i2==0;
eq9=delta(2)*i2+theta(2)*s2-gamma*v2==0;
eq10=i2+v2+s2+e2==1;

eq11=gamma*v3-theta(3)*s3-s3*(bE(3)*e2+bI(3)*i2)==0;
eq12=s3*(bE(3)*e2+bI(3)*i2)-epsilon*e3==0;
eq13=epsilon*e3-delta(3)*i3==0;
eq14=delta(3)*i3+theta(3)*s3-gamma*v3==0;
eq15=i3+v3+s3+e3==1;

eq16=gamma*v4-theta(4)*s4-s4*(bE(4)*e3+bI(4)*i3)==0;
eq17=s4*(bE(4)*e3+bI(4)*i3)-epsilon*e4==0;
eq18=epsilon*e4-delta(4)*i4==0;
eq19=delta(4)*i4+theta(4)*s4-gamma*v4==0;
eq20=i4+v4+s4+e4==1;

eq21=gamma*v1-theta(1)*s1-s1*(bE(5)*e1+bI(5)*i1+bE(5)*e3+bI(5)*i3)==0;
eq22=s5*(bE(5)*e1+bI(5)*i1+bE(5)*e3+bI(5)*i3)-epsilon*e5==0;
eq23=epsilon*e5-delta(5)*i5==0;
eq24=delta(5)*i5+theta(5)*s5-gamma*v5==0;
eq25=i5+v5+s5+e5==1;

[pv1,pv2,pv3,pv4,pv5,pe1,pe2,pe3,pe4,pe5,pi1,pi2,pi3,pi4,pi5,ps1,ps2,ps3,ps4,ps5]=solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq15,eq16,eq17,eq18,eq19,eq20,eq21,eq22,eq23,eq24,eq25,v1,v2,v3,v4,v5,e1,e2,e3,e4,e5,i1,i2,i3,i4,i5,s1,s2,s3,s4,s5);
epi_result=double([pv1,pv2,pv3,pv4,pv5,pe1,pe2,pe3,pe4,pe5,pi1,pi2,pi3,pi4,pi5,ps1,ps2,ps3,ps4,ps5]);
epi_result=epi_result';
P_star=1/N*(pe1+pe2+pe3+pe4+pe5+pi1+pi2+pi3+pi4+pi5)








% T=1000;
% v1=zeros(T,1);
% v2=zeros(T,1);
% v3=zeros(T,1);
% v4=zeros(T,1);
% v5=zeros(T,1);
% dv1=zeros(T,1);
% dv2=zeros(T,1);
% dv3=zeros(T,1);
% dv4=zeros(T,1);
% dv5=zeros(T,1);
% 
% s1=zeros(T,1);
% s2=zeros(T,1);
% s3=zeros(T,1);
% s4=zeros(T,1);
% s5=zeros(T,1);
% ds1=zeros(T,1);
% ds2=zeros(T,1);
% ds3=zeros(T,1);
% ds4=zeros(T,1);
% ds5=zeros(T,1);
% 
% e1=zeros(T,1);
% e2=zeros(T,1);
% e3=zeros(T,1);
% e4=zeros(T,1);
% e5=zeros(T,1);
% de1=zeros(T,1);
% de2=zeros(T,1);
% de3=zeros(T,1);
% de4=zeros(T,1);
% de5=zeros(T,1);
% 
% i1=zeros(T,1);
% i2=zeros(T,1);
% i3=zeros(T,1);
% i4=zeros(T,1);
% i5=zeros(T,1);
% di1=zeros(T,1);
% di2=zeros(T,1);
% di3=zeros(T,1);
% di4=zeros(T,1);
% di5=zeros(T,1);
% 
% s1(1)=0.9;
% s2(1)=0.9;
% s3(1)=0.9;
% s4(1)=0.9;
% s5(1)=0.9;
% 
% i1(1)=0.1;
% i2(1)=0.1;
% i3(1)=0.1;
% i4(1)=0.1;
% i5(1)=0.1;
% 
% for tt=1:T
% ds1(tt)=gamma*v1(tt)-theta(1)*s1(tt)-s1(tt)*(bE(1)*e5(tt)+bI(1)*i5(tt));
% de1(tt)=s1(tt)*(bE(1)*e5(tt)+bI(1)*i5(tt))-epsilon*e1(tt);
% di1(tt)=epsilon*e1(tt)-delta(1)*i1(tt);
% dv1(tt)=delta(1)*i1(tt)+theta(1)*s1(tt)-gamma*v1(tt);
% 
% ds2(tt)=gamma*v2(tt)-theta(2)*s2(tt)-s2(tt)*(bE(2)*e1(tt)+bI(2)*i1(tt)+bE(2)*e3(tt)+bI(2)*i3(tt)+bE(2)*e4(tt)+bI(2)*i4(tt));
% de2(tt)=s2(tt)*(bE(2)*e1(tt)+bI(2)*i1(tt)+bE(2)*e3(tt)+bI(2)*i3(tt)+bE(2)*e4(tt)+bI(2)*i4(tt))-epsilon*e2(tt);
% di2(tt)=epsilon*e2(tt)-delta(2)*i2(tt);
% dv2(tt)=delta(2)*i2(tt)+theta(2)*s2(tt)-gamma*v2(tt);
% 
% 
% ds3(tt)=gamma*v3(tt)-theta(3)*s3(tt)-s3(tt)*(bE(3)*e2(tt)+bI(3)*i2(tt));
% de3(tt)=s3(tt)*(bE(3)*e2(tt)+bI(3)*i2(tt))-epsilon*e3(tt);
% di3(tt)=epsilon*e3(tt)-delta(3)*i3(tt);
% dv3(tt)=delta(3)*i3(tt)+theta(3)*s3(tt)-gamma*v3(tt);
% 
% 
% ds4(tt)=gamma*v4(tt)-theta(4)*s4(tt)-s4(tt)*(bE(4)*e3(tt)+bI(4)*i3(tt));
% de4(tt)=s4(tt)*(bE(4)*e3(tt)+bI(4)*i3(tt))-epsilon*e4(tt);
% di4(tt)=epsilon*e4(tt)-delta(4)*i4(tt);
% dv4(tt)=delta(4)*i4(tt)+theta(4)*s4(tt)-gamma*v4(tt);
% 
% 
% ds5(tt)=gamma*v1(tt)-theta(1)*s1(tt)-s1(tt)*(bE(5)*e1(tt)+bI(5)*i1(tt)+bE(5)*e3(tt)+bI(5)*i3(tt));
% de5(tt)=s5(tt)*(bE(5)*e1(tt)+bI(5)*i1(tt)+bE(5)*e3(tt)+bI(5)*i3(tt))-epsilon*e5(tt);
% di5(tt)=epsilon*e5(tt)-delta(5)*i5(tt);
% dv5(tt)=delta(5)*i5(tt)+theta(5)*s5(tt)-gamma*v5(tt);
% 
% v1(tt+1)=v1(tt)+dv1(tt);
% v2(tt+1)=v2(tt)+dv2(tt);
% v3(tt+1)=v3(tt)+dv3(tt);
% v4(tt+1)=v4(tt)+dv4(tt);
% v5(tt+1)=v5(tt)+dv5(tt);
% 
% i1(tt+1)=i1(tt)+di1(tt);
% i2(tt+1)=i2(tt)+di2(tt);
% i3(tt+1)=i3(tt)+di3(tt);
% i4(tt+1)=i4(tt)+di4(tt);
% i5(tt+1)=i5(tt)+di5(tt);
% 
% e1(tt+1)=e1(tt)+de1(tt);
% e2(tt+1)=e2(tt)+de2(tt);
% e3(tt+1)=e3(tt)+de3(tt);
% e4(tt+1)=e4(tt)+de4(tt);
% e5(tt+1)=e5(tt)+de5(tt);
% 
% s1(tt+1)=s1(tt)+ds1(tt);
% s2(tt+1)=s2(tt)+ds2(tt);
% s3(tt+1)=s3(tt)+ds3(tt);
% s4(tt+1)=s4(tt)+ds4(tt);
% s5(tt+1)=s5(tt)+ds5(tt);
% 
% 
% 
% 
% 
% 
% end

sbI=sum(bI)
sbE=sum(bE)
sdelta=sum(delta)
stheta=sum(theta)








function hI_b = hI(beta)
UbetaI=0.6;
hI_b=1/beta-1/UbetaI;
end

function hE_b = hE(beta)
UbetaE=0.7;
hE_b=1/beta-1/UbetaE;
end

function gi = g(phi,delta)
Udelta=0.1;
gi=inv_pos(phi-delta)-1/(phi-Udelta);
end

function gih = ghat(phi,delta)
Udelta=0.1;
gih=inv_pos(delta)-1/(phi-Udelta);
end

function fi = f(gamma,theta)
Utheta=0.1;
fi=(theta+gamma)/gamma-(Utheta+gamma)/gamma;
end

function fih = fhat(gamma,tao)
Utheta=0.1;
fih=inv_pos(tao)-1-Utheta/gamma;
end



