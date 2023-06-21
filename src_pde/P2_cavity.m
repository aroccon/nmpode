%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

%-----------------------------------------------------------------------
Re = 2e2;     % Reynolds number
dt = 1e-3;    % time step
tf = 10e-0;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
nx = 90;      % number of x-gridpoints
ny = 90;      % number of y-gridpoints
nsteps = 10;  % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;           % number of time steps
x = linspace(0,lx,nx+1); hx = lx/nx;    % Grid along x
y = linspace(0,ly,ny+1); hy = ly/ny;    % Grid along y
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx-1,ny); V = zeros(nx,ny-1);
% boundary conditions (N=North-top, 
uN = x*0+1;    vN = avg(x)*0;
uS = x*0;      vS = avg(x)*0;
uW = avg(y)*0; vW = y*0;
uE = avg(y)*0; vE = y*0;
%-----------------------------------------------------------------------
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
      [uW;zeros(nx-3,ny);uE]/hy^2);
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);

fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:nt
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   % Extension on boundary point
   Ue = [uW;U;uE]; 
   Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; 
   Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   %x-derivative of UV
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   %y-derivative of UV
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg(Ue(:,2:end-1));   
   Ud = diff(Ue(:,2:end-1))/2;
   Va = avg(Ve(2:end-1,:)')'; 
   Vd = diff(Ve(2:end-1,:)')'/2;
   
   %computer d(u^2)/dx and d(v^2)/dy 
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   % Update values of interior points:
   U = U-dt*(UVy(2:end-1,:)+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);
   rhs = reshape(V+Vbc,[],1);
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   % Compute divergence of the velocity field
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   % Update values of interior points by 
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy; 
end

figure(1)
clf
subplot(1,2,1)
hold on
imagesc(x,y,Ue')
axis equal
title('X - Velocity')
hold off


subplot(1,2,2)
hold on
title('Y - Velocity')
imagesc(x,y,Ve')
axis equal
hold off




function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end
end

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
end

