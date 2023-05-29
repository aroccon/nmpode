dt = [ 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5]; 
t_final = 1;
for j = 1:length(dt)
  [t,y] =  P1_midpointsolver(1,dt(j),t_final,@P1_derivs1st); 
  X(j,2) = abs(y(end,1) - exp(-1));
  X(j,1) = dt(j);
end
loglog(X(:,1),X(:,2))
hold on 
loglog(X(:,1),X(:,1).^2,'--r')