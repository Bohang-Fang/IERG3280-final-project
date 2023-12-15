% Model parameters
beta = 5*10^-9; % rate of infection
gamma = 0.12; % rate of recovery (try also 0.07)
delta = 0; % rate of immunity loss
imloss=60; %immunity is lost after 60 days
N = 6*10^7; % Total population N = S + I + R
I0 = 10; % initial number of infected
T = 300; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)
fprintf('Value of parameter R0 is %.2f',(N-I0)*beta/gamma)
% Calculate the model
[S,I,R] = sirs_model(beta,gamma,delta,N,I0,T,dt);
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
figure(1)
plot(tt,S,'b','LineWidth',2);hold on
plot(tt,I,'r','LineWidth',2);hold on
plot(tt,R,'g','LineWidth',2);hold on;grid on;
xlabel('Days','FontSize',14); ylabel('Number of individuals','FontSize',14);
title('SIRS model');
legend('S','I','R');
legend('FontSize',12)
hold off
% Map
% figure(2)
% plot(I(1:(T/dt)-1),I(2:T/dt),"LineWidth",1,"Color",'r');
% hold on; grid on;
% plot(I(2),I(1),'ob','MarkerSize',4);
% xlabel('Infected at time t'); ylabel('Infected at time t+1');
% hold off;

[S1,I1,R1] = sirl_model(beta,gamma,imloss,N,I0,T,dt);
figure(2)
plot(tt,S1,'b','LineWidth',2);hold on
plot(tt,I1,'r','LineWidth',2);hold on
plot(tt,R1,'g','LineWidth',2);hold on;grid on;
xlabel('Days','FontSize',14); ylabel('Number of individuals','FontSize',14);
title('Recovered state all changes to susceptible state');
legend('S','I','R');
legend('FontSize',12)




function [S,I,R] = sirs_model(beta,gamma,delta,N,I0,T,dt)
    % if delta = 0 we have a SIR model (particular case)
    S = zeros(1,T/dt);
    S(1) = N-I0;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    for tt = 1:(T/dt)-1
        % Equations of the model
        dS = (-beta*I(tt)*S(tt) + delta*R(tt)) * dt;
        dI = (beta*I(tt)*S(tt) - gamma*I(tt)) * dt;
        dR = (gamma*I(tt) - delta*R(tt)) * dt;
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
    end
end

function [S,I,R] = sirl_model(beta,gamma,imloss,N,I0,T,dt)
    % if delta = 0 we have a SIR model (particular case)
    S = zeros(1,T/dt);
    S(1) = N-I0;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    for tt = 1:(T/dt)-1
        if tt<=imloss/dt-1;
        % Equations of the model
            dS(tt) = (-beta*I(tt)*S(tt) ) * dt;
            dI(tt) = (beta*I(tt)*S(tt) - gamma*I(tt)) * dt;
            dR(tt) = (gamma*I(tt) ) * dt; 
            S(tt+1) = S(tt) + dS(tt);
            I(tt+1) = I(tt) + dI(tt);
            R(tt+1) = R(tt) + dR(tt);
        else tt>imloss/dt-1;
            dS(tt) = (-beta*I(tt)*S(tt) ) * dt;
            dI(tt) = (beta*I(tt)*S(tt) - gamma*I(tt)) * dt;
            dR(tt) = (gamma*I(tt) ) * dt;
            S(tt+1) = S(tt) + dS(tt)+dR(tt-imloss/dt+1);
            I(tt+1) = I(tt) + dI(tt);
            R(tt+1) = R(tt) + dR(tt)-dR(tt-imloss/dt+1);
        end
        
    end
end