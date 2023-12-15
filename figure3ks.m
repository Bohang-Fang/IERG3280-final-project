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
k=0;

cvx_begin gp
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


    minimize(obj)
    subject to
        
        epsilon*u(1)+deltah(1)*v(1)<=(phi-k)*v(1);
        epsilon*u(2)+deltah(2)*v(2)<=(phi-k)*v(2);
        epsilon*u(3)+deltah(3)*v(3)<=(phi-k)*v(3);
        epsilon*u(4)+deltah(4)*v(4)<=(phi-k)*v(4);
        epsilon*u(5)+deltah(5)*v(5)<=(phi-k)*v(5);

        tao(1)*(bE(1)*u(5)+bI(1)*v(5))+(phi-epsilon)*u(1)<=(phi-k)*u(1);
        tao(2)*(bE(2)*u(1)+bI(2)*v(1))+tao(2)*(bE(2)*u(3)+bI(2)*v(3))+tao(2)*(bE(2)*u(4)+bI(2)*v(4))+(phi-epsilon)*u(2)<=(phi-k)*u(2);
        tao(3)*(bE(3)*u(2)+bI(3)*v(2))+(phi-epsilon)*u(3)<=(phi-k)*u(3);
        tao(4)*(bE(4)*u(4)+bI(4)*v(4))+(phi-epsilon)*u(4)<=(phi-k)*u(4);
        tao(5)*(bE(5)*u(1)+bI(5)*v(1))+tao(5)*(bE(5)*u(3)+bI(5)*v(3))+(phi-epsilon)*u(5)<=(phi-k)*u(5);
        
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


obj/C

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