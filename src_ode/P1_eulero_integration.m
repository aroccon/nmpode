% Explicit Eulero solver
clear all

% time step
dx=0.1;
x_init=0;
x_final=2;
Nsteps = round(x_final/dx); 

% define x array
x = linspace(0,x_final,Nsteps+1);
xth = linspace(0,x_final,100);
yth = exp(-xth);
y = zeros(1,length(x));

% Initial condition
y(1) = 1.;


for i = 1:Nsteps
      y(i+1) = y(i) - dx*y(i);
end


figure(1)
clf
hold on
plot(xth,yth,'--r','DisplayName','Analytical solution');
set(gca,'FontSize',30)
plot(x,y,'Displayname','Numerical solution');
legend show
hold off


