function [p] = P2_poissonfast(nx,rhs)

% 3D FFT for Fourier-spectral solution of
% p_xx + p_yy  + p_zz = rhs;

% wavenumbers
k = [0:(nx/2-1) (-nx/2):(-1)]; 
[kx,ky,kz]  = meshgrid(k,k,k);
% Laplacian matrix acting on the wavenumbers
delsq = -(kx.^2 + ky.^2 + kz.^2); 
% Laplacian matrix acting on the wavenumbers
% avodi solving zero waevnumber
delsq(1,1,1) = 1;  

% trasnform of rhsp
fhat = fftn(rhs);
p = real(ifftn(fhat./delsq));
% Specify arbitrary constant by forcing corner p = 0.
p = p - p(1,1,1);   


end