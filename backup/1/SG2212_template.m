% Navier-Stokes solver,
% adapted for course SG2212
% KTH Mechanics
%
% Depends on avg.m and DD.m
%
% Code version:
% 20180222
 
clear all;
close all;
clc;

%----------------(A)---------------------

lid_driven_cavity=1;

if (lid_driven_cavity==1)
% Parameters for test case I: Lid-driven cavity
% The Richardson number is zero, i.e. passive scalar.

Pr = 0.71; % Prandtl number
Re = 4000; % Reynolds number
Ri = 0.; % Richardson number

dt = 0.001; % time step
Tf = 20; % final time
Lx = 1; % width of box
Ly = 1; % height of box
Nx = 100; % number of cells in x
Ny = 100; % number of cells in y
ig = 200; % number of iterations between output

% Boundary and initial conditions:
Utop = 1.;
Ubottom = 0.;
% IF TEMPERATURE: Tbottom = 1.; Ttop = 0.;
% IF TEMPERATURE: namp = 0.;
else
% Parameters for test case II: Rayleigh-Brd convection
% The DNS will be stable for Ra=1705, and unstable for Ra=1715
% (Theoretical limit for pure sinusoidal waves
% with L=2.01h: Ra=1708)
% Note the alternative scaling for convection problems.

Pr = 0.71; % Prandtl number
Ra = 1708; % Rayleigh number

Re = 1./Pr; % Reynolds number
Ri = Ra*Pr; % Richardson number

dt = 0.001; % time step
Tf = 20; % final time
Lx = 10.; % width of box
Ly = 1; % height of box
Nx = 30; % number of cells in x
Ny = 30; % number of cells in y
ig = 200; % number of iterations between output

% Boundary and initial conditions:
Utop = 0.;
Ubottom = 0.;
% IF TEMPERATURE: Tbottom = 1.; Ttop = 0.;
namp = 0.1;
end


%-----------------------------------------

% Number of iterations
Nit = Tf/dt;
% Spatial grid: Location of corners
x = linspace( 0,Lx,Nx+1 );
y = linspace( 0,Ly,Ny+1 );
% Grid spacing

dx = Lx/Nx;
dy = Ly/Ny;
% Boundary conditions:
uN = x*0+Utop; vN = avg(x,2)*0;
uS = zeros(1,Nx+1); vS = zeros(1,Nx);
uW = zeros(1,Ny); vW = zeros(1,Ny+1);
uE = zeros(1,Ny); vE = zeros(1,Ny+1);
tN = linspace(0,Tf,Tf/dt); tS = linspace(0,Tf,Tf/dt);
% Initial conditions
U = zeros(Nx-1,Ny);  V = zeros(Nx,Ny-1);
% linear profile for T with random noise
% IF TEMPERATURE: T = ... + namp*rand(Nx,Ny)
% Time series
tser = [];
Tser = [];

%-----------------------------------------

% Compute system matrices for pressure
% First set homogeneous Neumann condition all around
% Laplace operator on cell centres: Fxx + Fyy
Lp = kron(speye(Ny), DD(Nx,dx) ) + kron( DD(Ny,dy) ,speye(Nx));
% Set one Dirichlet value to fix pressure in that point
Lp(1,:) = 0 ; Lp(1,1)= 1 ;
% Here you can pre-compute the LU decomposition
% [LLp,ULp] = lu(Lp);
%-----------------------------------------

% Progress bar (do not replace the ... )
fprintf(...
'[           |           |           |           |           ]\n')

%-----------------------------------------

% Main loop over iterations

for k = 1:Nit

% include all boundary points for u and v (linear extrapolation
% for ghost cells) into extended array (Ue,Ve)
Ue = [uW;U;uE];
Ue= [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
Ve = [2*vS'-V(:,1) V 2*vN'-V(:,end)];
Ve=[vW; Ve; vE];

% averaged (Ua,Va) of u and v on corners
Ua = avg(Ue,2);
Va = avg(Ve,1);

% construct individual parts of nonlinear terms
dUVdx = diff(Ua.*Va,1,1)/dx;
dUVdy = diff(Ua.*Va,1,2)/dy;
dU2dx = diff(Ua.*Ua,1,1)/dx;
dV2dy = diff(Va.*Va,1,2)/dy;

% treat viscosity explicitly
viscu = diff(Ue(:,2:end-1),2,1)/dx^2+diff(Ue(2:end-1,:),2,2)/dy^2;
viscv = diff(Ve(:,2:end-1),2,1)/dx^2+diff(Ve(2:end-1,:),2,2)/dy^2;

% buoyancy term
% IF TEMPERATURE: fy = ...

% compose final nonlinear term + explicit viscous terms
U = U + dt/Re*viscu - dt*(dU2dx(1:end-1,2:end)+dUVdy(2:end-1,:));
V = V + dt/Re*viscv - dt*(dUVdx(:,2:end-1)+dV2dy(1:end-1,2:end)); % IF TEMPERATURE: dt*fy;

% pressure correction, Dirichlet P=0 at (1,1)
rhs = (diff( [uW;U;uE] )/dx + diff([2*vS'-V(:,1) V 2*vN'-V(:,end)],1,2)/dy)/dt;
rhs = reshape(rhs,Nx*Ny,1);
rhs(1) = 0;
P = Lp\rhs;
% alternatively, you can use the pre-computed LU decompositon
% P = ...;
% or gmres
% P = gmres(Lp, rhs, [], tol, maxit);
% or as another alternative you can use GS / SOR from homework 6
% [PP, r] = GS_SOR(omega, Nx, Ny, hx, hy, L, f, p0, tol, maxit);
P = reshape(P,Nx,Ny);

% apply pressure correction
U = U - dt*diff(P,1,1)/dx;
V = V - dt*diff(P,1,2)/dy;

% Temperature equation
% IF TEMPERATURE: Te = ...
% IF TEMPERATURE: Tu = ...
% IF TEMPERATURE: Tv = ...
% IF TEMPERATURE: H = ...
% IF TEMPERATURE: T = T + dt*H;

%-----------------------------------------

% progress bar
if floor(51*k/Nit)>floor(51*(k-1)/Nit), fprintf('.'), end

% plot solution if needed
if k==1|floor(k/ig)==k/ig

% compute divergence on cell centres
if (1==1)
div = diff([uW;U;uE])/dx + diff([vS' V vN'],1,2)/dy;

figure(1);clf; hold on;
contourf(avg(x,2),avg(y,2),div');colorbar
axis equal; axis([0 Lx 0 Ly]);
title(sprintf('divergence at t=%g',k*dt))
drawnow
end

% compute velocity on cell corners
Ua = avg(Ue,2);
Va = avg(Ve,1);
Len = sqrt(Ua.^2+Va.^2+eps);

figure(2);clf;hold on;
%contourf(avg(x,2),avg(y,2),P');colorbar
contourf(x,y,sqrt(Ua.^2+Va.^2)',20,'k-');colorbar
quiver(x,y,(Ua./Len)',(Va./Len)',.4,'k-')
axis equal; axis([0 Lx 0 Ly]);
title(sprintf('u at t=%g',k*dt))
drawnow

% IF TEMPERATURE: % compute temperature on cell corners
% IF TEMPERATURE: Ta = ...

% IF TEMPERATURE: figure(3); clf; hold on;
% IF TEMPERATURE: contourf(x,y,Ta',20,'k-');colorbar
% IF TEMPERATURE: quiver(x,y,(Ua./Len)',(Va./Len)',.4,'k-')
% IF TEMPERATURE: axis equal; axis([0 Lx 0 Ly]);
% IF TEMPERATURE: title(sprintf('T at t=%g',k*dt))
% IF TEMPERATURE: drawnow

% Time history
if (1==1)
figure(4); hold on;
tser = [tser k*dt];
Tser = [Tser Ue(ceil((Nx+1)/2),ceil((Ny+1)/2))];
plot(tser,abs(Tser))
title(sprintf('Probe signal at x=%g, y=%g',...
x(ceil((Nx+1)/2)),y(ceil((Ny+1)/2))))
set(gca,'yscale','log')
xlabel('time t');ylabel('u(t)')
end
end
end
fprintf('\n')