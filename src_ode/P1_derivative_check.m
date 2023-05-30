clear all
% define time array
x = linspace(0,2,20);
y = exp(-x);

% numerical derivative
dydx = diff(y)./diff(x);

figure(1)
clf
hold on
plot(x,y,'DisplayName','f(x)=exp(x)');
set(gca,'FontSize',30)
plot(x(1:end-1),dydx,'r--','DisplayName','Numeric derivative');
plot(x,-y,'Displayname','Analytical derivative');
legend show
hold off


