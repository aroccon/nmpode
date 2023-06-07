% Gauss-Seidel iterative solver for Poisson equation
% 07/07/2023 A. Roccon
% Remark: code written assuming dx=dy, bear this in mind if changing lx or
% ly or the grid resolution.
clear all

%input parameters
nx=32;
ny=nx;
lx=1.0;
ly=1.0;
max_err=1.e-4;
max_iter=40000;
err=2.0d0;
dx=lx/(nx-1);
dy=dx;
x=zeros(nx,1);
y=zeros(ny,1);



for i=1:nx-1
    x(i+1)=x(i)+dx;
end
for j=1:ny-1
    y(j+1)=y(j)+dy;
end
    

% start
A=zeros(nx,ny);
A0=zeros(nx,ny);
An=zeros(nx,ny);
rhs=zeros(nx,ny);

%Boundary conditions (temperature gradient on left and bottom sides).
for i=1:nx
  A(i,1)=1.;%double(i-1)/double(nx-1);
end
for j=2:ny
  A(nx,j)=0.;%double(ny-j)/double(ny-1);
end

% Store initial condition (just to check)
A0=A;
% Updtae An with Boundary conditons, otherwise new values are not updated.
An=A0;

% Computer RHS
% Generic function
 for i=2:nx-1
   for j=2:ny-1
      %rhs(i,j)=10*(x(i)^2 + y(j)^2 + 10);
      rhs(i,j)=10*(1-x(i)^2);
   end
end
 
%source in the center of the domain
%rhs(nx/2,ny/2)=-400;

iter=1;
while ((err >= max_err) && ( iter <= max_iter))
  for i=2:nx-1
    for j=2:ny-1
       An(i,j)=0.25*(A(i+1,j)+An(i-1,j)+A(i,j+1)+An(i,j-1)-dx*dx*rhs(i,j));
    end
  end
  % Compute error
  err=0.0d0;
  for i=2:nx-1
     for j=2:ny-1
       err = max(abs(An(i,j) - A(i,j)),err);
       A(i,j) = An(i,j);
     end
  end
  iter = iter+1;
end

figure(1)
clf
subplot(1,2,1)
hold on
imagesc(x,flipud(y),A0);
title('Initial condition')
hold off


subplot(1,2,2)
hold on
title('Final solution')
imagesc(x,flipud(y),A);
hold off

