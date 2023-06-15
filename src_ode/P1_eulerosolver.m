%% Euler Solver function
function [t,data] = P1_eulersolver(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt);
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; 
data(1,:) = y';
for i =1:Nsteps
     dy = feval(derivs_Handle,time,y);
     y = y + dt*dy;
     time = time + dt;
     t(i+1) = time;
     data(i+1,:)  = y';
     plot(time,y(1),'*');
end

dy = feval(derivs_Handle,time,y); 
y = y + dy*dt;
time = time+dt;
t(i+1) = time;
data(i+1,:)  = y';