% Jacobi solver for Laplace equation
% 07/07/2023 A. Roccon
% No assumption is made on dx=dy, more flexible code.
clear all

%input parameters
nx=32;
ny=128;
lx=1.0;
ly=2.0;
max_err=1.e-4;
max_iter=10000;
err=2.0d0;
dx=lx/(nx-1);
dy=ly/(ny-1);
x=zeros(nx,1);
y=zeros(ny,1);
coeff=(2./dx/dx + 2./dy/dy);



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

%Boundary conditions (temperature gradient on left and bottom sides).
for i=1:nx
  A(i,1)=double(i-1)/double(nx-1);
end
for j=1:ny
  A(nx,j)=double(ny-j)/double(ny-1);
end

% Store initial condition (just to check)
A0=A;

iter=1;
while ((err >= max_err) && ( iter <= max_iter))
  for i=2:nx-1
    for j=2:ny-1
       An(i,j)=1./(coeff*dx*dx)*(A(i+1,j)+A(i-1,j)) + 1./(coeff*dy*dy)*(A(i,j+1)+A(i,j-1));
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
subplot(1,2,1)
hold on
imagesc(x,flipud(y),A0)
axis equal
title('Initial condition')
hold off


subplot(1,2,2)
hold on
title('Final solution')
imagesc(x,flipud(y),A)
axis equal
hold off

