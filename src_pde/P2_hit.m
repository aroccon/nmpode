% A. Roccon 08/02/2024
% Homogenous isotropic turbulence solver
% Constant density and viscosity
% 2nd order finite difference + fastPoisson solver for pressure 
% ABC forcing scheme (see Comparison of forcing schemes to sustain
% homogeneous isotropic turbulence)
% Must be validated

clear
% parameters
% domain size is 2*pi along all directions
% same grid spacing
lx=2*pi;
nx=32; % for pressure without boundary nodes
dx=2*pi/nx; 
tfin=20000;
dt=0.001;
dxi=1/dx;
ddxi=1/dx/dx;
rho=1; %rho
nu=0.005; %mu/rho

% forcing parameters (usually A=B=C)
A=0.2;
B=0.2;
C=0.2;
k0=5;

%define axis of staggered grid (staggered)
xs=linspace(-lx/2,lx/2,nx);
x=linspace(-lx/2,lx/2,nx+1); %only inner points

% Variable defined on cell centers (pressure)
rhsp=zeros(nx,nx,nx);
p=zeros(nx,nx,nx);

% Variables defined on cell faces
u=zeros(nx+1,nx+1,nx+1);
v=zeros(nx+1,nx+1,nx+1);
w=zeros(nx+1,nx+1,nx+1);
div=zeros(nx+1,nx+1,nx+1);
fx=zeros(nx+1,nx+1,nx+1);
fy=zeros(nx+1,nx+1,nx+1);
fz=zeros(nx+1,nx+1,nx+1);
ustar=zeros(nx+1,nx+1,nx+1);
vstar=zeros(nx+1,nx+1,nx+1);
wstar=zeros(nx+1,nx+1,nx+1);
rhsu=zeros(nx+1,nx+1,nx+1);
rhsv=zeros(nx+1,nx+1,nx+1);
rhsw=zeros(nx+1,nx+1,nx+1);



% Init velocity field
u(:,:,:)=0;
v(:,:,:)=0;
w(:,:,:)=0;

disp('Starting temporal loop')

% Start temporal loop
for t=1:tfin
disp(dt*t)

% Projection step, convective terms


% Convective terms NS
% Compute here the products h11=u1*u1 h12=u1*u2 h21=u2*u1 and h22=u2*u2
% Average between the values at i+1*i and i-1*i
for i=1:nx+1
  for j=1:nx+1
      for k=1:nx+1
        % periodiciy without requiring extranodes
        ip=i+1;
        im=i-1;
        jp=j+1;
        jm=j-1;
        kp=k+1;
        km=k-1;
        if (ip > nx); ip=1; end
        if (im < 1);  im=nx; end
        if (jp > nx); jp=1; end
        if (jm < 1);  jm=nx; end        
        if (kp > nx); kp=1; end
        if (km < 1);  km=nx; end
        % compute the products (conservative form)
        h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k));
        h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jp,k))*(v(i,j,k)+v(im,j,k));
        h13 = (u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k));
        h21 = (u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k));
        h22 = (v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k));
        h23 = (w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km));
        h31 = (w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km));
        h32 = (v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k));
        h33 = (w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km));
        % compute the derivative
        h11=h11*0.25*dxi;
        h12=h12*0.25*dxi;
        h13=h13*0.25*dxi;
        h21=h21*0.25*dxi;
        h22=h22*0.25*dxi;
        h23=h23*0.25*dxi;
        h31=h31*0.25*dxi;
        h32=h32*0.25*dxi;
        h33=h33*0.25*dxi;
        % add to the rhs
        rhsu(i,j,k)=-(h11+h12+h13);
        rhsv(i,j,k)=-(h21+h22+h23);
        rhsw(i,j,k)=-(h31+h32+h33);
      end
  end
end

% Compute viscous terms
for i=1:nx+1
    for j=1:nx+1
        for k=1:nx+1
            % periodiciy without requiring extranodes
            ip=i+1;
            im=i-1;
            jp=j+1;
            jm=j-1;
            kp=k+1;
            km=k-1;
            if (ip > nx); ip=1; end
            if (im < 1);  im=nx; end
            if (jp > nx); jp=1; end
            if (jm < 1);  jm=nx; end        
            if (kp > nx); kp=1; end
            if (km < 1);  km=nx; end
            % compute diffusive terms (second order derivative)
            h11 = nu*(u(ip,j,k)-2*u(i,j,k)+u(im,j,k))*ddxi;
            h12 = nu*(u(i,jp,k)-2*u(i,j,k)+u(i,jm,k))*ddxi;
            h13 = nu*(u(i,j,kp)-2*u(i,j,k)+u(i,j,km))*ddxi;
            h21 = nu*(v(ip,j,k)-2*v(i,j,k)+v(im,j,k))*ddxi;
            h22 = nu*(v(i,jp,k)-2*v(i,j,k)+v(i,jm,k))*ddxi;
            h23 = nu*(v(i,j,kp)-2*v(i,j,k)+v(i,j,km))*ddxi;
            h31 = nu*(w(ip,j,k)-2*w(i,j,k)+w(im,j,k))*ddxi;
            h32 = nu*(w(i,jp,k)-2*w(i,j,k)+w(i,jm,k))*ddxi;
            h33 = nu*(w(i,j,kp)-2*w(i,j,k)+w(i,j,km))*ddxi;
            rhsu(i,j,k)= rhsu(i,j,k)+(h11+h12+h13);
            rhsv(i,j,k)= rhsv(i,j,k)+(h21+h22+h23);
            rhsw(i,j,k)= rhsw(i,j,k)+(h31+h32+h33);
        end
    end
end

%forcing term (always x because is the same axis)
for i=1:nx+1
    for j=1:nx+1
        for k=1:nx+1
            fx(i,j,k)=A*sin(k0*x(k))+C*sin(k0*x(j));
            fy(i,j,k)=B*sin(k0*x(i))+A*sin(k0*x(k));
            fz(i,j,k)=C*sin(k0*x(j))+B*sin(k0*x(i));
            rhsu(i,j,k)= rhsu(i,j,k) + fx(i,j,k);
            rhsv(i,j,k)= rhsv(i,j,k) + fy(i,j,k);
            rhsw(i,j,k)= rhsw(i,j,k) + fz(i,j,k);
        end
    end
end

% find u, v and w star (explicit Eulero)
for i=1:nx+1
    for j=1:nx+1
        for k=1:nx+1
            ustar(i,j,k) = u(i,j,k) + dt*rhsu(i,j,k);
            vstar(i,j,k) = v(i,j,k) + dt*rhsv(i,j,k);
            wstar(i,j,k) = w(i,j,k) + dt*rhsw(i,j,k);
        end
    end
end


% Compute rhs of Poisson equation div*ustar
% Compute divergence at the cell center (nx+1 avilable on u,v,w)
for i=1:nx
    for j=1:nx
        for k=1:nx
            rhsp(i,j,k) = (rho*dxi/dt)*(ustar(i+1,j,k)-ustar(i,j,k) + vstar(i,j+1,k)-vstar(i,j,k) + wstar(i,j,k+1)-wstar(i,j,k));
        end
    end
end

% call Poisson solver (3DFastPoissonsolver, periodic BCs)
% Input only boundary nodes
p = P2_poissonfast(nx,rhsp);


% Correct velocity 
for i=1:nx
    for j=1:nx
        for k=1:nx
            % periodicity 
            im=i-1;
            jm=j-1;
            km=k-1;
            if (im < 1);  im=nx; end
            if (jm < 1);  jm=nx; end        
            if (km < 1);  km=nx; end
            u(i,j,k)=ustar(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi;
            v(i,j,k)=vstar(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dxi;
            w(i,j,k)=wstar(i,j,k) - dt/rho*(p(i,j,k)-p(i,j,km))*dxi;
        end
    end
end

% Check divergence
for i=1:nx
    for j=1:nx
        for k=1:nx
            div(i,j,k) = (u(i+1,j,k)-u(i,j,k))*dxi + (v(i,j+1,k)-v(i,j,k))*dxi + (w(i,j,k+1)-w(i,j,k))*dxi;
        end
    end
end

% update halo and impose periodicity
% x-periodicity 
u(1,:,:)=u(nx+1,:,:);
v(1,:,:)=v(nx+1,:,:);
w(1,:,:)=w(nx+1,:,:);
% y-periodicity 
u(:,1,:)=u(:,nx+1,:);
v(:,1,:)=v(:,nx+1,:);
w(:,1,:)=w(:,nx+1,:);
% z-periodicity 
u(:,:,1)=u(:,:,nx+1);
v(:,:,1)=v(:,:,nx+1);
w(:,:,1)=w(:,:,nx+1);


% end of temporal loop
end


vtkwrite('u.vtk', 'structured_points', 'u',u(2:end-1,2:end-1,2:end-1))
%vtkwrite('v.vtk', 'structured_points', 'v',v(2:end-1,2:end-1,2:end-1))
%vtkwrite('w.vtk', 'structured_points', 'w',w(2:end-1,2:end-1,2:end-1))
%vtkwrite('p.vtk', 'structured_points', 'p',p(2:end-1,2:end-1,2:end-1))
