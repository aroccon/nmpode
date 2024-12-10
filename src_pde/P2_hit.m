% A. Roccon 08/02/2024
% Homogenous isotropic turbulence solver
% Constant density and viscosity
% 2nd order finite difference + fastPoisson solver for pressure 
% ABC forcing scheme (see Comparison of forcing schemes to sustain
% homogeneous isotropic turbulence)
% Qualitative results seem fine

clear all
% parameters
% domain size is 2*pi along all directions
% same grid spacing
lx=2*pi;
nx=64; % for pressure without boundary nodes
dx=2*pi/(nx-1); 
tfin=50000;
dt=0.001;
dxi=1/dx;
ddxi=1/dx/dx;
rho=1; %rho
mu=0.03; %mu/rho

% forcing parameters (Lundgren)
A=1;
B=1;
C=1;
k0=2;


%define axis of staggered grid (staggered)
x=linspace(-lx/2,lx/2,nx);
%x=linspace(-lx/2,lx/2,nx+1); %only inner points

% Variable defined on cell centers (pressure)
rhsp=zeros(nx,nx,nx);
p=zeros(nx,nx,nx);

% Variables defined on cell faces
u=zeros(nx,nx,nx);
v=zeros(nx,nx,nx);
w=zeros(nx,nx,nx);
div=zeros(nx,nx,nx);
fx=zeros(nx,nx,nx);
fy=zeros(nx,nx,nx);
fz=zeros(nx,nx,nx);
ustar=zeros(nx,nx,nx);
vstar=zeros(nx,nx,nx);
wstar=zeros(nx,nx,nx);
rhsu=zeros(nx,nx,nx);
rhsv=zeros(nx,nx,nx);
rhsw=zeros(nx,nx,nx);


% Init velocity field
% Taylor-Green + Random noise
%u=2*(u-0.5);
%v=2*(v-0.5);
%w=2*(w-0.5);

for k = 1:nx
    for j= 1:nx
        for i = 1:nx
            u(i,j,k) =  sin(x(i)) * cos(x(j)) * cos(x(k));
            v(i,j,k) =  cos(x(i)) * sin(x(j)) * cos(x(k));
            w(i,j,k) =  0;
        end
    end
end

disp('Starting temporal loop')

% Start temporal loop
for t=1:tfin
%disp(dt*t)

% Projection step, convective terms


% Convective terms NS
for i=1:nx
  for j=1:nx
      for k=1:nx
        ip=i+1;
        jp=j+1;
        kp=k+1;
        im=i-1;
        jm=j-1;
        km=k-1;
        if (ip > nx); ip=1; end
        if (jp > nx); jp=1; end
        if (kp > nx); kp=1; end   
        if (im < 1); im=nx; end
        if (jm < 1); jm=nx; end
        if (km < 1); km=nx; end   
        % compute the products (conservative form)
        h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k));
        h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k));
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
for i=1:nx
    for j=1:nx
        for k=1:nx
            % compute diffusive terms (second order derivative)
            ip=i+1;
            jp=j+1;
            kp=k+1;
            im=i-1;
            jm=j-1;
            km=k-1;
            if (ip > nx); ip=1; end
            if (jp > nx); jp=1; end
            if (kp > nx); kp=1; end
            if (im < 1); im=nx; end
            if (jm < 1); jm=nx; end
            if (km < 1); km=nx; end   
            h11 = mu*(u(ip,j,k)-2*u(i,j,k)+u(im,j,k))*ddxi;
            h12 = mu*(u(i,jp,k)-2*u(i,j,k)+u(i,jm,k))*ddxi;
            h13 = mu*(u(i,j,kp)-2*u(i,j,k)+u(i,j,km))*ddxi;
            h21 = mu*(v(ip,j,k)-2*v(i,j,k)+v(im,j,k))*ddxi;
            h22 = mu*(v(i,jp,k)-2*v(i,j,k)+v(i,jm,k))*ddxi;
            h23 = mu*(v(i,j,kp)-2*v(i,j,k)+v(i,j,km))*ddxi;
            h31 = mu*(w(ip,j,k)-2*w(i,j,k)+w(im,j,k))*ddxi;
            h32 = mu*(w(i,jp,k)-2*w(i,j,k)+w(i,jm,k))*ddxi;
            h33 = mu*(w(i,j,kp)-2*w(i,j,k)+w(i,j,km))*ddxi;
            rhsu(i,j,k)= rhsu(i,j,k)+(h11+h12+h13)/rho;
            rhsv(i,j,k)= rhsv(i,j,k)+(h21+h22+h23)/rho;
            rhsw(i,j,k)= rhsw(i,j,k)+(h31+h32+h33)/rho;
        end
    end
end


%forcing term (always x because is the same axis)
for i=1:nx
    for j=1:nx
        for k=1:nx
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
for i=1:nx
    for j=1:nx
        for k=1:nx
            ustar(i,j,k) = u(i,j,k) + dt*rhsu(i,j,k);
            vstar(i,j,k) = v(i,j,k) + dt*rhsv(i,j,k);
            wstar(i,j,k) = w(i,j,k) + dt*rhsw(i,j,k);
        end
    end
end

%disp(max(max(max(ustar))))
%disp(max(max(max(vstar))))
%disp(max(max(max(wstar))))

% Compute rhs of Poisson equation div*ustar
% Compute divergence at the cell center (nx+1 avilable on u,v,w)
for i=1:nx
    for j=1:nx
        for k=1:nx
            ip=i+1;
            jp=j+1;
            kp=k+1;
            if (ip > nx); ip=1; end
            if (jp > nx); jp=1; end
            if (kp > nx); kp=1; end   
            rhsp(i,j,k) = (rho*dxi/dt)*(ustar(ip,j,k)-ustar(i,j,k));
            rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(vstar(i,jp,k)-vstar(i,j,k));
            rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(wstar(i,j,kp)-wstar(i,j,k));
        end
    end
end

%disp(max(max(max(rhsp))))

% call Poisson solver (3DFastPoissonsolver, periodic BCs)
p = P2_fastPoisson3D(nx,rhsp);


% Correct velocity 
for i=1:nx
    for j=1:nx
        for k=1:nx
            im=i-1;
            jm=j-1;
            km=k-1;
            if (im < 1); im=nx; end
            if (jm < 1); jm=nx; end
            if (km < 1); km=nx; end   
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
            ip=i+1;
            jp=j+1;
            kp=k+1;
            if (ip > nx); ip=1; end
            if (jp > nx); jp=1; end
            if (kp > nx); kp=1; end   
            div(i,j,k) = dxi*(u(ip,j,k)-u(i,j,k) + v(i,jp,k)-v(i,j,k) + w(i,j,kp)-w(i,j,k));
        end
    end
end

disp(max(max(max(div))))

%check courant number
uc=max(max(max(u)));
vc=max(max(max(v)));
wc=max(max(max(w)));

cou=max(uc*dt*dxi,vc*dt*dxi);
cou=max(cou,wc*dt*dxi);
if isnan(cou)
    disp('NaN..!!')
    break
end



% end of temporal loop
end


vtkwrite('u.vtk', 'structured_points', 'u',u)
vtkwrite('v.vtk', 'structured_points', 'v',v)
vtkwrite('w.vtk', 'structured_points', 'w',w)
vtkwrite('p.vtk', 'structured_points', 'p',p)