% 2D FFT for Fourier-spectral solution of
% u_xx + u_yy  = (r^2 - 2*sig^2)*exp(-r^2/(2*sig^2))/sig^4,  0 < x,y < 1
% where
% r^2 = (x-0.5)^2 + (y-0.5)^2, sig = 0.2
% uex = exp(-r^2/(2*sig^2)) + C (C is an arbitrary constant
% We will choose C by taking u(x=0,y=0) = 0. In the limit sig << 1, this
% implies C = 0.

% No. of Fourier modes...should be a power of 2
N = 32;      
% Domain size (assumed square
L = 1;      
sig = 0.1;   

% wavenumbers
k = (2*pi/L)*[0:(N/2-1) (-N/2):(-1)]; 
[KX, KY]  = meshgrid(k,k);
% Laplacian matrix acting on the wavenumbers
delsq = -(KX.^2 + KY.^2); % Laplacian matrix acting on the wavenumbers
% avodi solving zero waevnumber
delsq(1,1) = 1;          

% Construct RHS f(x,y) at the Fourier gridpoints
h = L/N;     % Grid spacing
x = (0:(N-1))*h ;
y = (0:(N-1))*h;
[X Y] = meshgrid(x,y);
rsq = (X-0.5*L).^2 + (Y-0.5*L).^2;
sigsq = sig^2;
f = exp(-rsq/(2*sigsq)).*(rsq - 2*sigsq)/(sigsq^2);

% trasnform of forcing
fhat = fft2(f);
u = real(ifft2(fhat./delsq));
% Specify arbitrary constant by forcing corner u = 0.
u = u - u(1,1);   

%  Plot out solution in interior

uex = exp(-rsq/(2*sigsq));
%errmax = norm(u(:)-uex(:),inf)
contourf(X,Y,u)
xlabel('x')
ylabel('y')
zlabel('u')
title('Fourier spectral method for 2D Poisson Eqn')