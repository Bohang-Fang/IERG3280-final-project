figure(1)
s=[0.5;0.6;0.7;0.8;0.9;1.0;1.1;1.2;1.3;1.4;1.5;1.75;2];
lamda_star=[0.139012;0.106547;0.0767911;0.0493078;0.0238325;0; -0.022421;-0.043549;-0.0634864;-0.0823261;-0.100169;-0.140965;-0.177062];
plot(s,lamda_star,'-r*','Linewidth',1.2,'MarkerSize',6);grid on;
xlabel('s','FontSize',14)
y=ylabel('$\lambda^*$','interpreter','latex','FontSize',14)
set(y,'Rotation',0)


figure(2)
k=[0;0.05;0.01;0.03;0.07;0.1;0.12;0.15];
S=[1;1.2317;1.0439;1.1352;1.3339;1.4990;1.6176;1.8098];
plot(k,S,'-r*','Linewidth',1.2,'MarkerSize',7);grid on;
xlabel('k','FontSize',14)

y2=ylabel('$\frac{obj}{C^*}$','interpreter','latex','FontSize',20)
axis([0,0.15,0.8,2])
set(y2,'Rotation',0)


figure(3)
s1=[0;0.1;0.2;0.3;0.7;1;1.3;1.7;2;2.5;3;4;5.5];
bI=[3;3;3;3;3;2.9186;2.7292; 2.4753;2.3038;2.0633;1.8665;1.5641;1.3434];
bE=[3.5;3.5;3.5;3.5;3.4671;3.3054;3.1637;2.9813;2.8378;2.6257;2.4414;2.1357;1.8920];
delta=[0.5;1.0597;1.2381;1.4495;2.0707;2.3372;2.5316;2.7375;2.8657;3.0450;3.1910;3.4137;3.5747];
theta=[0.5;0.519;0.5871;0.6494;0.9;1.0583;1.1905;1.3561; 1.4741;1.6725;1.8730;2.2801;2.6953];

plot(s1,bI,'-sb','Linewidth',1.2,'MarkerSize',6); hold on;grid on
plot(s1,bE,'-or','Linewidth',1.2,'MarkerSize',6); hold on;
plot(s1,delta,'-g+','Linewidth',1.2,'MarkerSize',6); hold on;
plot(s1,theta,'-k.','Linewidth',1.2,'MarkerSize',12); hold on;
xlabel('s','FontSize',20)
axis([0,5.5,0.49,4])

legend({'$\sum \beta^I$','$\sum \beta^E$','$\sum \delta$','$\sum \theta$'},'interpreter','latex')
legend('FontSize',11)
set(gca,'FontName','Times New Roman')






