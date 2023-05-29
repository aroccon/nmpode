clear all;
[T,Y] = P1_midpointsolver([1;0],0.1,20,@derivs); 

figure(1)
clf 
subplot(1,2,1)
hold on
title('Position')
set(gca,'FontSize',30)
plot(T,cos(T),'r','linewidth',2.0,'Displayname','Analytical'); 
plot(T,Y(:,1),'--b','linewidth',2.0,'Displayname','Numerical');  
legend show
hold off

subplot(1,2,2)
hold on
title('Velocity')
set(gca,'FontSize',30)
plot(T,-sin(T),'r','linewidth',2.0,'Displayname','Analytical'); 
plot(T,Y(:,2),'--b','linewidth',2.0,'Displayname','Numerical'); 
legend show
hold off