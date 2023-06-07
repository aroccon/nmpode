%Implementation of a vorticity-streamfunction method for quasi 2D-turbulence.
%Space discretization: Pseudo-spectral.
%Time-discretization: Third-order Runge-Kutta 

M = 256; % number of points
N = M;
Lx = 2*pi;
Ly = 2*pi;
nu = 5e-4; % kinematic viscosity 
Sc = 0.7; % Schmidt number
ar = 0.12; %random number amplitude
b = 1; % mean scalar gradient
CFLmax = 0.8; 
tend = 200; % end time

x=linspace(0,Lx,M+1); x(end)=[];
dx = Lx/M;
kx=[0:M/2 -M/2+1:-1]*2*pi/Lx;

y=linspace(0,Ly,N+1); y(end)=[];
ky=[0:N/2 -N/2+1:-1]*2*pi/Ly;
dy = Ly/N;

time = 0;

index_kmax = ceil(M/3);
kmax = kx(index_kmax);
filter = ones(M,N);
filter(index_kmax+1:2*index_kmax+3,index_kmax+1:2*index_kmax+3)=0;

rng(64);

[u, v, omega, psi, ddx, ddy, idel2, kk, k2]=deal(zeros(M,N));

%Kx-Wavenumber
for j=1:N
    ddx(:,j)=1i*kx;
end

%Ky_wavenumber
for i=1:M
    ddy(i,:)=1i*ky;
end

% Wavenumbers squared 
for i=1:M
    for j=1:N
        idel2(i,j)=-kx(i)^2-ky(j)^2;
    end
end
% Inverse of wavenumbers squared 1/|k^2|
idel2=1./idel2;
idel2(1,1)=0;

for i=1:M
    for j=1:N
        kk(i,j)=kx(i)^2+ky(j)^2;
        k2(i,j)=kx(i)^2+ky(j)^2;
        if kk(i,j) >= 6^2 && kk(i,j) <= 7^2 % forcing 
           kk(i,j) = -kk(i,j); 
        end
        
        if kk(i,j) <= 2^2
            kk(i,j) = 8*kk(i,j); % large-scale dissipation
        end
    end
end

% Initialize velocity: Taylor-Green + Random noise
for i=1:M
    for j=1:N
        u(i,j) =  cos(2*x(i))*sin(2*y(j))+ar*rand;
        v(i,j) = -sin(2*x(i))*cos(2*y(j))+ar*rand;
    end
end

% u and v in spectral space
uc = fft2(u);
vc = fft2(v);
% Compute omega (vorticity) in spectral space 
omegac = ddx.*vc - ddy.*uc; 

% Initialize phi with random noise + bring it to the spectral space% (phic)
phi = rand(size(u));
phic = fft2(phi);

dt = 0.5*min([dx dy]);

nstep = 1;

while time < tend
    % Compute streamfunction \nabla^2 \psi = -omega and from psi, uc and vc
    psic = -idel2.*omegac;     
    uc = ddy.*psic;           
    vc = -ddx.*psic;           
    
    % Compute u, v from uc, vc
    u = real(ifft2(uc));       
    v = real(ifft2(vc));      
    
    % Compute the derivates of omega w.r.t. to x and y.
    omegadx = real(ifft2(ddx.*omegac));
    omegady = real(ifft2(ddy.*omegac));
    
    facto = exp(-nu*8/15*dt*kk);
    factp = exp(-nu/Sc*8/15*dt*k2);
    
    r0o = -fft2(u.*omegadx+v.*omegady);
    r0p = -fft2(u.*real(ifft2(ddx.*phic))+v.*real(ifft2(ddy.*phic)))+b*vc;
    
    omegac = facto.*(omegac + dt*8/15*r0o); % update omega
    phic = factp.*(phic + dt*8/15*r0p); % update phi

    %%%% Substep 2
    psic = -idel2.*omegac;
    uc = ddy.*psic;
    vc = -ddx.*psic;
    
    u = real(ifft2(uc));
    v = real(ifft2(vc));
    
    omegadx = real(ifft2(ddx.*omegac));
    omegady = real(ifft2(ddy.*omegac));
    
    r1o = -fft2(u.*omegadx+v.*omegady);
    r1p = -fft2(u.*real(ifft2(ddx.*phic))+v.*real(ifft2(ddy.*phic)))+b*vc;
    
    omegac = omegac + dt*(-17/60*facto.*r0o + 5/12*r1o);
    phic = phic + dt*(-17/60*factp.*r0p + 5/12*r1p);
    facto = exp(-nu*(-17/60+5/12)*dt*kk);
    factp = exp(-nu/Sc*(-17/60+5/12)*dt*k2);
    omegac = omegac.*facto;
    phic = phic.*factp;
    
    %%%% Substep 3
    psihat = -idel2.*omegac;
    uc = ddy.*psihat;
    vc = -ddx.*psihat;
    
    % max(max(abs(real(ifft2(1i*ddx.*uc+1i*ddy.*vc))))) % divergence
    
    u = real(ifft2(uc));
    v = real(ifft2(vc));
    
    omegadx = real(ifft2(ddx.*omegac));
    omegady = real(ifft2(ddy.*omegac));
    
    r2o = -fft2(u.*omegadx+v.*omegady);
    r2p = -fft2(u.*real(ifft2(ddx.*phic))+v.*real(ifft2(ddy.*phic)))+b*vc;    
    omegac = omegac + dt*(-5/12*facto.*r1o + 3/4*r2o);
    phic = phic + dt*(-5/12*factp.*r1p + 3/4*r2p);
    facto = exp(-nu*(-5/12+3/4)*dt*kk);
    factp = exp(-nu/Sc*(-5/12+3/4)*dt*kk);
    omegac = omegac.*facto;
    phic = phic.*factp;

    phic = filter.*phic;
    omegac = filter.*omegac;
    
    time = time + dt;
    nstep = nstep + 1;
    
    CFL = max(max(abs(u)))/dx*dt+max(max(abs(v)))/dy*dt;
    
    if mod(nstep,20)==0
        phi = real(ifft2(phic));
        omega = real(ifft2(omegac));
        dissipation = 2*nu*(real(ifft2(ddx.*uc)).^2 + real(ifft2(ddy.*uc)).^2 + real(ifft2(ddx.*vc)).^2 + real(ifft2(ddy.*vc)).^2);
        eta = (nu^3/mean(dissipation,'all'))^0.25;
        
        subplot(221); pcolor(x,y,omega'); title('Vorticity'); shading flat; axis equal tight; colorbar; drawnow
        subplot(222); pcolor(x,y,phi'); title('Scalar');  shading flat; axis equal tight; colorbar; drawnow
        subplot(223); pcolor(x,y,dissipation'); title('Dissipation'); shading flat; axis equal tight; colorbar; drawnow
        subplot(224); pcolor(x,y,u); title('u-velocity'); shading flat; axis equal tight; colorbar; drawnow
        
        fprintf(1,'step = %d    time = %g    dt = %g  CFL = %g    kmax*eta = %g %g\n', nstep, time, dt, CFL, eta.*kmax, eta.*kmax/sqrt(Sc));
        
    end
    
    dt = CFLmax/CFL*dt; %0.0005;
end

