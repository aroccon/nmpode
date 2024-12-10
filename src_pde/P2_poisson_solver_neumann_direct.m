function [p] = P2_poisson_solver_neumann_direct(xs,ys,rhs)

% The location of p（NOTE: Staggered Grid）
% Input from the calling subroutine
dx=xs(2)-xs(1);
dy=ys(2)-ys(1);
%[x,y] = ndgrid(xs,ys); % Note: not meshgrid
nx=length(xs);
ny=length(ys);

% Constructing a matrix A (sparse matrix to save memory usage)
% Start with the x-derivative (tri-diagonal matrix)
tmp = -2*ones(nx,1);
tmp([1,end]) = -1; % Neumann BC requires -1 (not -2) at both ends
Ad = diag(tmp);
Au = diag(ones(nx-1,1),1);
Al = diag(ones(nx-1,1),-1);
Ax = Ad+Au+Al;

% Consider y-derivative = it become block diagnal matrix
% Define the container first
% Abig = zeros(Nx*Ny,Nx*Ny,'like',sparse(1));
% dd = eye(Nx);
% for jj=1:Ny
%     if jj==1 || jj==Ny % y-ends（Neumann BC)
%         Abig(1+Nx*(jj-1):Nx*jj,1+Nx*(jj-1):Nx*jj) = Ax/dx^2 - dd/dy^2;
%     else
%         Abig(1+Nx*(jj-1):Nx*jj,1+Nx*(jj-1):Nx*jj) = Ax/dx^2 - 2*dd/dy^2;
%     end
%     if jj~=1 % j
%         Abig(1+Nx*(jj-1):Nx*jj,1+Nx*(jj-2):Nx*(jj-1)) = dd/dy^2;
%     end
%     if jj~=Ny
%         Abig(1+Nx*(jj-1):Nx*jj,1+Nx*(jj-0):Nx*(jj+1)) = dd/dy^2;
%     end
% end

% The above can be written a little simpler:
% Construct the block diagonal matrix
dd = eye(nx);
tmp = repmat({sparse(Ax/dx^2 - 2*dd/dy^2)},ny,1);
tmp{1} = Ax/dx^2 - dd/dy^2;
tmp{ny} = Ax/dx^2 - dd/dy^2; % y-ends（Neumann BC)
Abig = blkdiag(tmp{:}); 

% A matrix (one block smaller)
d4y = eye(nx*(ny-1),'like',sparse(1));
Abig(1:end-nx,nx+1:end) = Abig(1:end-nx,nx+1:end) + d4y/dy^2; % upper
Abig(nx+1:end,1:end-nx) = Abig(nx+1:end,1:end-nx) + d4y/dy^2; % lower

% rhs is from input
f=reshape(rhs,[],1);

% Natually, Abig is a singular matrix. Poisson eq with Neumann BC
% does not provide a unique solution. Thus fixing u to be 0 at one point.
Abig(1,:) = 0;
Abig(1,1) = 1;
f(1) = 0;

% Solve the system of equation
p = Abig\f;

% Put the solution (a vector) back to 2d matrix
p = reshape(p,[nx,ny]);
%p = fliplr(p');
%p = fliplr(p);

end