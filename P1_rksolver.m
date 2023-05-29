%% Runge-Kutta Solver
function [t,data] = P1_rksolver(y,dt,t_final,P1_derivs_Handle) 
time = 0;
Nsteps = round(t_final/dt);  
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition data(1,:) = y’;
for i =1:Nsteps
    k1 = dt*feval(P1_derivs_Handle,time ,y ); 
    k2 = dt*feval(P1_derivs_Handle,time+dt/2,y+k1/2);
    k3 = dt*feval(P1_derivs_Handle,time+dt/2,y+k2/2); 
    k4=dt*feval(P1_derivs_Handle,time+dt ,y+k3 ); 
    y = y + k1/6 + k2/3 + k3/3 + k4/6;
    time = time+dt;
    t(i+1) = time;
    data(i+1,:)  = y';
end