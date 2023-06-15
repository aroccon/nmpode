%% Implicit Eulero Solver function
function [t,data] = P1_impeulerosolver(y,dt,t_final)
time = 0;
Nsteps = round(t_final/dt);
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; 
data(1,:) = y';

% Exploit the fact that the mass-spring system is linear avoiding iterative methods
% Reformulate the problem in a matrix form
A=zeros(2,2);
A=[1, -dt; dt, 1];

for i =1:Nsteps
     y = inv(A)*y;
     time = time + dt;
     t(i+1) = time;
     data(i+1,:)  = y';
     plot(time,y(1),'*');
end
