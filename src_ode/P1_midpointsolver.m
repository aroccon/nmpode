function [t,data] = P1_midpointsolver(y,dt,t_final,derivs_Handle) 
time = 0;
Nsteps = round(t_final/dt);
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; 
data(1,:) = y';
for i =1:Nsteps
    %% evaluate the initial derivatives
    dy = feval(derivs_Handle,time,y);  
    yH = y + dy*dt/2; %% take Euler step to midpoint
    dy = feval(derivs_Handle,time,yH); 
    y = y + dy*dt; 
    time = time+dt; 
    t(i+1) = time; 
    data(i+1,:) = y';
end