% Gauss-Seidel iterative solver for Poisson equation
% 07/07/2023 A. Roccon
% Remark: code written assuming dx=dy, bear this in mind if changing lx or
% ly or the grid resolution.
% Pressure field for a static droplet
clear all

%input parameters
nx=64;
ny=nx;
lx=2.0;
ly=2.0;
max_err=1.e-4;
max_iter=40000;
ch=0.04; % Transition layer thickness for diffuse interface droplet
radius=0.5; % Drop radius
sigma=1.0;

err=2.0d0;
dx=lx/(nx-1);
dy=dx;
x=zeros(nx,1)-lx/2;
y=zeros(ny,1)-ly/2;

for i=1:nx-1
    x(i+1)=x(i)+dx;
end
for j=1:ny-1
    y(j+1)=y(j)+dy;
end
    
% start + boundary conditions on pressure (0 at the boundaries)
A=zeros(nx,ny);
A0=zeros(nx,ny);
An=zeros(nx,ny);

% Computer RHS (surface tension forces)
% Interfacial layer thickness
phi=zeros(nx,ny);
for i=1:nx
    for j=1:ny
        phi(i,j)=-tanh((sqrt(x(i)^2 +y(j)^2)-radius)/sqrt(2)/ch);
    end
end
        
gphix=zeros(nx,ny);
gphiy=zeros(nx,ny);
mgphi=zeros(nx,ny);

% compute gradient and mangnitude
for i=2:nx-1
    for j=2:ny-1
        gphix(i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dx);
        gphiy(i,j)=(phi(i,j+1)-phi(i,j-1))/(2*dy);
        mgphi(i,j)=sqrt(gphix(i,j)^2+gphiy(i,j)^2)+1.e-9; 
    end
end

% compute curvature
k=zeros(nx,ny);
normx=zeros(nx,ny);
normy=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        normx(i,j)=gphix(i,j)/mgphi(i,j);
        normy(i,j)=gphiy(i,j)/mgphi(i,j);
    end
end
for i=2:nx-1
    for j=2:ny-1
        k(i,j)= - (normx(i+1,j)-normx(i-1,j))/(2*dx) - (normy(i,j+1)-normy(i,j-1))/(2*dy);
    end
end

% compute surface tension forces
fx=zeros(nx,ny);
fy=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        fx(i,j)=3/sqrt(8)*ch*sigma*k(i,j)*gphix(i,j)*mgphi(i,j);
        fy(i,j)=3/sqrt(8)*ch*sigma*k(i,j)*gphiy(i,j)*mgphi(i,j);
    end
end

% computer rhs (divergence of fx, fy)
rhs=zeros(nx,ny);
for i=5:nx-4
    for j=5:ny-4
        rhs(i,j)= - (fx(i+1,j)-fx(i-1,j))/(2*dx) - (fy(i,j+1)-fy(i,j-1))/(2*dy);
    end
end



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
subplot(2,2,1)
hold on
imagesc(x,flipud(y),phi);
title('Droplet shape')
hold off

subplot(2,2,2)
hold on
imagesc(x,flipud(y),rhs);
title('Right-hand side')
hold off

subplot(2,2,3)
hold on
imagesc(x,flipud(y),A0);
title('Initial condition')
hold off


subplot(2,2,4)
hold on
title('Final solution')
imagesc(x,flipud(y),A);
hold off

