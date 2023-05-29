% First test Euler explicit for mass spring
% Not optimized from a computational point of view (figure always open)
clear;
y(1) = 1; 
y(2) = 0;
dt = 0.1;
time = 0;
t_final = 30; 
Nsteps = round(t_final/dt);
plot(time,y(1),'*');
hold on

for i = 1:Nsteps
    dy(2)  = -y(1);
    dy(1)  =  y(2);
    y = y + dt*dy;
    time = time + dt;
    plot(time,y(1),'b*');
    plot(time,cos(time),'r.');
end