% A. Roccon 08/02/2024
% Validated with re=100 + dt=0.001
% Constant density and viscosity 
% Primitive variables and projection-correction
% Poisson solver with pure Neumann conditions (dp/dn=0 everywhere is solid boundary)
% Based on direct solver (Abig inversion)
% Pressure imposed on one point to avoid singular matrix

clear all
% parameters
lx=1.0;
ly=1.0;
nx=64;
ny=64;
tfin=5000;
dt=0.001;
dx=lx/(nx);
dy=ly/(ny);
dxi=1/dx;
dyi=1/dy;
dti=1/dt;
re=100;
mu=1.0;
rho=1;


%define axis of staggered grid (staggered)
xs=linspace(-lx/2,lx/2,nx+2);
ys=linspace(-ly/2,ly/2,ny+2);
x=linspace(-lx/2,lx/2,nx+3);
y=linspace(-ly/2,ly/2,ny+3);

% Variable defined on cell centers (phi and pressure)
phi=zeros(nx+2,ny+2);
p=zeros(nx+2,ny+2);
rhsp=zeros(nx+2,ny+2);


% Variables defined on cell faces
u=zeros(nx+3,ny+3);
v=zeros(nx+3,ny+3);
ustar=zeros(nx+3,ny+3);
vstar=zeros(nx+3,ny+3);
tau12=zeros(nx+3,ny+3);
h11=zeros(nx+3,ny+3);
h12=zeros(nx+3,ny+3);
h21=zeros(nx+3,ny+3);
h22=zeros(nx+3,ny+3);



% Init velocity field
u(:,:)=0.0;
v(:,:)=0.0;

% 
u(:,ny+2:ny+3)=1.0;



% Init phase-field
for i=1:nx+1
    for j=1:ny+1
        pos=xs(i)^2+ys(j)^2;
        phi(i,j)=0;%0.5*(1-tanh((sqrt(pos)-radius)/2/eps));
    end
end

disp('Staring temporal loop')
disp(eps/dx)


% Start temporal loop
for t=1:tfin
    disp(t)

% Projection step, convective terms
% Compute rhs 
rhsu=zeros(nx+3,ny+3);
rhsv=zeros(nx+3,ny+3);


% Convective terms NS
% Compute here the products h11=u1*u1 h12=u1*u2 h21=u2*u1 and h22=u2*u2
% Average between the values at i+1*i and i-1*i
for i=2:nx+2
  for j=2:ny+2
      % compute the products (conservative form)
      h11 = (u(i+1,j)+u(i,j))*(u(i+1,j)+u(i,j))- (u(i,j)+u(i-1,j))*(u(i,j)+u(i-1,j));
      h12 = (u(i,j+1)+u(i,j))*(v(i,j+1)+v(i-1,j+1))- (u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j));
      h21 = (u(i+1,j)+u(i+1,j-1))*(v(i+1,j)+v(i,j))- (u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j));
      h22 = (v(i,j+1)+v(i,j))*(v(i,j+1)+v(i,j))-(v(i,j)+v(i,j-1))*(v(i,j)+v(i,j-1));
      % compute the derivative
      h11=h11*0.25*dxi;
      h12=h12*0.25*dyi;
      h21=h21*0.25*dxi;
      h22=h22*0.25*dyi;
      % add to the rhs
      rhsu(i,j)=h11+h12;
      rhsv(i,j)=h21+h22;
  end
end

% Compute viscous terms
% Assume uniform viscosity and density
for i=2:nx+3
    for j=2:ny+3
        tau12(i,j) = mu*((u(i,j)-u(i,j-1))*dyi+(v(i,j)-v(i-1,j))*dxi);
    end
end

for i=2:nx+2
    for j=2:ny+2
         % compute diffusive terms
        h11 = (1/re)*(u(i+1,j)-2*u(i,j)+u(i-1,j))*dxi*dxi;
        h22 = (1/re)*(v(i,j+1)-2*v(i,j)+v(i,j-1))*dyi*dyi;
        h12 = (1/re)*(tau12(i,j+1)-tau12(i,j))*dyi;
        h21 = (1/re)*(tau12(i+1,j)-tau12(i,j))*dxi;
        % equivalent to the second order derivative i+1 -2i +i-1
        rhsu(i,j)= rhsu(i,j)-(h11+h12);
        rhsv(i,j)= rhsv(i,j)-(h22+h21);
    end
end

% find u and v star (explicit Eulero)
for i=2:nx+2
    for j=2:ny+2
        ustar(i,j) = u(i,j) - dt*rhsu(i,j);
        vstar(i,j) = v(i,j) - dt*rhsv(i,j);
    end
end


% impose BC on the star field
ustar(1:2,:)=0.0;
ustar(nx+2:nx+3,:)=0.0;
ustar(:,1:2)=0.0;
ustar(:,ny+2:ny+3)=1.0;

vstar(:,1:2)=0.0;
vstar(:,nx+2:nx+3)=0.0;
vstar(1:2,:)=0.0;
vstar(ny+2:ny+3,:)=0.0;



% Compute rhs of Poisson equation div*ustar
for i=1:nx+2
    for j=1:ny+2
     rhsp(i,j) = (rho/dt)*((ustar(i+1,j)-ustar(i,j))*dxi + (vstar(i,j+1)-vstar(i,j))*dyi);
    end
end

% call Poisson solver
% It is a mess to use density variatons and Neumann BCs, because it is not
% so easy to assembel the matrix Abig (see function poisson solver).

p=poisson_solver_neumann_direct(xs,ys,rhsp);

% Correct velocity 
for i=2:nx+2
    for j=2:ny+2
        u(i,j)=ustar(i,j) - dt/rho*(p(i,j)-p(i-1,j))*dxi;
        v(i,j)=vstar(i,j) - dt/rho*(p(i,j)-p(i,j-1))*dyi;
    end
end

% Check divergence
for i=1:nx+2
    for j=1:ny+2
     div(i,j) = (u(i+1,j)-u(i,j))*dxi + (v(i,j+1)-v(i,j))*dyi;
    end
end

% impose BC on the new field
u(1:2,:)=0.0;
u(nx+2:nx+3,:)=0.0;
u(:,1:2)=0.0;
u(:,ny+2:ny+3)=1.0;

v(:,1:2)=0.0;
v(:,nx+2:nx+3)=0.0;
v(1:2,:)=0.0;
v(ny+2:ny+3,:)=0.0;


% end of temporal loop
end


figure(1)
subplot(2,3,1)
contourf(x(2:end-1),y(2:end-1),u(2:end-1,2:end-1)',20)
hold on
colorbar
title('U Velocity')
hold off

subplot(2,3,2)
contourf(x(2:end-1),y(2:end-1),v(2:end-1,2:end-1)',20)
hold on
colorbar
title('V Velocity')
hold off

subplot(2,3,3)
contourf(xs(2:end-1),ys(2:end-1),p(2:end-1,2:end-1)',20)
hold on
colorbar
title('Pressure')
hold off

subplot(2,3,4)
contourf(xs(2:end-1),ys(2:end-1),div(2:end-1,2:end-1)',20)
hold on
title('Divergence (check)')
colorbar
hold off

subplot(2,3,5)
contourf(x(2:end-1),y(2:end-1),(sqrt(u(2:end-1,2:end-1).^2 + v(2:end-1,2:end-1).^2))',20)
hold on
title('Velocity Magnitude')
colormap('turbo')
colorbar
hold off

