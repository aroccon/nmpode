
clear all

% ti


%error of subroutines P1 computed as
% yth(final) - ynum(final)
dx=[0.001,0.005,0.01,0.05,0.1,0.5];
err=[1.35e-4,6.77e-4,0.0014,0.0058,0.0138,0.103];
scaling=2*dx.^1
scaling2=2*dx.^2


figure(1)
clf
loglog(dx,err,'--r','DisplayName','Error');
hold on
loglog(dx,scaling,'b','DisplayName','Order P=1');
loglog(dx,scaling2,'k','DisplayName','Order P=2');
hold on
set(gca,'FontSize',30)
legend show
hold off



