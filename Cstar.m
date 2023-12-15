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



cvx_begin gp
    variables bI(5) 
    variables bE(5) 
    variables delta(5)
    variables tao(5) 
    variables u(5)
    obj=0;

    for ii=1:N
         obj=obj+hI(bI(ii))+hE(bE(ii))+g(phi,delta(ii))+fhat(gamma,tao(ii));
    end



    minimize(obj)
    subject to
        tao(1)*u(5)*(bE(1)+epsilon/delta(5)*bI(1))+(phi-epsilon)*u(1)<=u(1)*(phi-barK);
        tao(2)*u(1)*(bE(2)+epsilon/delta(1)*bI(2))+tao(2)*u(3)*(bE(2)+epsilon/delta(3)*bI(2))+tao(2)*u(4)*(bE(2)+epsilon/delta(4)*bI(2))+(phi-epsilon)*u(2)<=u(2)*(phi-barK);
        tao(3)*u(2)*(bE(3)+epsilon/delta(2)*bI(3))+(phi-epsilon)*u(3)<=u(3)*(phi-barK);
        tao(4)*u(3)*(bE(4)+epsilon/delta(3)*bI(4))+(phi-epsilon)*u(4)<=u(4)*(phi-barK);
        tao(5)*u(1)*(bE(5)+epsilon/delta(1)*bI(5))+tao(5)*u(3)*(bE(5)+epsilon/delta(3)*bI(5))+(phi-epsilon)*u(5)<=u(5)*(phi-barK);
        gamma/(1+gamma)<=tao(1)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(2)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(3)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(4)<=gamma/(0.1+gamma);
        gamma/(1+gamma)<=tao(5)<=gamma/(0.1+gamma);
        0.1<=delta(1)<=1;
        0.1<=delta(2)<=1;
        0.1<=delta(3)<=1;
        0.1<=delta(4)<=1;
        0.1<=delta(5)<=1;
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


theta=gamma*(1-tao)./tao;
C_star=0;

for ii=1:N
    C_star=C_star+hI(bI(ii))+hE(bE(ii))+g(phi,delta(ii))+f(gamma,theta(ii));
end

BI=diag(bI);
BE=diag(bE);
E=epsilon*eye(N);
D=diag(delta);
T=diag(gamma*ones(N,1)./(theta+gamma*ones(N,1)));
Q=[T*BE*A-E T*BI*A;
    E -D];

max(real(eig(Q)))


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

function gih = gh(phi,delta)
Udelta=0.1;
gih=1/(delta)-1/(phi-Udelta);
end

function fi = f(gamma,theta)
Utheta=0.1;
fi=(theta+gamma)/gamma-(Utheta+gamma)/gamma;
end

function fih = fhat(gamma,tao)
Utheta=0.1;
fih=inv_pos(tao)-1-Utheta/gamma;
end