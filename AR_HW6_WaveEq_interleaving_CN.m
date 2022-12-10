%% MAE207 SIMULATION ODE's HW6 - Q1
%  Andres Rodriguez A53275066   12/08/2022
%  Solve Hyperbolic Wave equation with pade in space CN in time
%  Similar solution as algorithm 11.6 but: 
%  Interleave the x=[q v] vector to march CN w/o iterations
%  Note: interleaving => x = [q1 v1 q2 v2 .... and so on ]
%  ====================================================================

close all; 
clear
% Definde discretization values for time and space
P=1; %1=> Periodic BC. 0=>Dirichlect BC
L = 4;
N=64;
dt = 0.01;
dx=L/N;
c=1;
Tmax = 4;

% Create matrices to be used for interleaving
I = eye(N);
La = 0.1*ones(N-1,1); Da = ones(N,1); Ua = La;
Lb = 1.2*ones(N-1,1); Db = -2.4*ones(N,1); Ub = Lb;
A = diag(La,-1)+diag(Da,0)+diag(Ua,1);
B = (diag(Lb,-1)+diag(Db,0)+diag(Ub,1))*(1/dx^2); 

clear La Lb Da Db Ua Ub
C1 = [I zeros(N,N);    zeros(N,N) A];
C2 = [zeros(N,N) I; (c^2)*B zeros(N,N)];


%% INTERLEAVE MATRICES
% Reorder C1,C2 => D1, D2
ctr1=1; ctr2=1; mid=N;
% Reorder rows
for i=1:length(C1)
    if rem(i,2)==1 %odd row
        temp1(i,:) = C1(ctr1,:);
        temp2(i,:) = C2(ctr1,:);
        ctr1=ctr1+1;
    elseif rem(i,2)==0 %even row
        temp1(i,:) = C1(ctr2+mid,:);
        temp2(i,:) = C2(ctr2+mid,:);
        ctr2=ctr2+1;
    end
end
% Reorder columns
ctr1=1; ctr2=1;
for i=1:length(temp1)
    if rem(i,2)==1 %odd row
        D1(:,i) = temp1(:,ctr1);
        D2(:,i) = temp2(:,ctr1);
        ctr1=ctr1+1;
    elseif rem(i,2)==0; %even row 
        D1(:,i) = temp1(:,ctr2+mid);
        D2(:,i) = temp2(:,ctr2+mid);
        ctr2=ctr2+1;
    end
end
%Adjust Corner Entries for Periodic Boundary condition
D1(2,2*N)=0.1; D1(2*N,2)=0.1; 
D2(2,2*N-1)=1.2/dx^2; D2(2*N,1)=1.2/dx^2;
C3 = inv(D1)*D2;

%% INITIALIZE ALGORTIHM
IterSteps=2; t=0; x=(-N/2:N/2-1)'*dx; q=exp(-x.^2/0.1); v=zeros(length(q),1);
PlotXY(x,q,t,-L/2,L/2,-0.2,1.2); vs=v; qs=v; a=0.6*dt*c^2/dx^2; b=-1.2*dt*c^2/dx^2;

%interleave q & v at t=0. => (y vector)
% Reorder rows
ctr1=1;ctr2=1;
for i=1:length(q)*2
    if rem(i,2)==1 %odd row
        y(i,1) = q(ctr1); ctr1=ctr1+1;
    elseif rem(i,2)==0 %even row
        y(i,:) = v(ctr2); ctr2=ctr2+1;
    end
end
% Extract diagonals from M3 H3 for the Penta Solver
% Solving only for interior gridpoints. Remove 1st and last 2 entries
L2=diag(D1,-2); L1=diag(D1,-1); D=diag(D1,0); U1=diag(D1,1); U2=diag(D1,2);

%% EXECUTE ALGORITHM 
% Solved with A\b rather than the Penta solver

if P==1 %(Periodic BC) -------------------------------------------------
    for n=1:Tmax/dt
    RHS = (1/dt*D1+0.5*D2)*y;
    yn1 = (1/dt*D1-0.5*D2)\RHS; %Ax=b found via A\b 
    %extract q values from yn+1 for plotting
    ctr=1; 
    for i=1:length(yn1)
       if rem(i,2)~= 0; qn1(ctr,1)=yn1(i); ctr=ctr+1; end
    end
    %update variables and plot
    t=t+dt; q=qn1; y=yn1; PlotXY(x,q,t,-L/2,L/2,-0.2,1.2);

    end
else %(Dirichlet BC) ---------------------------------------------------
    % Non-periodic BC. Algorithm Solves for interior points only 
    % Remove 2 top & Bottom rows of M3 H3 matrices
    D1(1:2,:)=[]; D1(end-1:end,:)=[]; D2(1:2,:)=[]; D2(end-1:end,:)=[];
    % Compute BC (since they depend on (t)
    %y(1) = y(end-1); %q(1) = q(end)
    %y(2) = y(end);   %v(1) = v(end)
    for n=1:Tmax/dt
    RHS = (1/dt*D1+0.5*D2)*y;
    yn1 = (1/dt*D1-0.5*D2)\RHS; %Ax=b found via A\b
    %Enforce BC
    %y(end-1)=0; y(end)=0;
    
    %extract q values from yn+1 for plotting
    ctr=1; 
    for i=1:length(yn1)
       if rem(i,2)~= 0; qn1(ctr,1)=yn1(i); ctr=ctr+1; end
    end
    %update variables and plot
    t=t+dt; q=qn1; y=yn1; PlotXY(x,q,t,-L/2,L/2,-1,1.2);
    end
end
set(gcf,'color','white')

%% PLOT XY

function PlotXY(x,y,t,xmin,xmax,ymin,ymax)
% function PlotXY(x,y,t,xmin,xmax,ymin,ymax)
% A supplemental plotting code used in several of the simulations in Chapter 11 of NR.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 11.2.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap11">Chapter 11</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.

plot(x,y);  xlabel('x');  ylabel('u'); grid on
title(sprintf('Time = %5.2f',t)); axis([xmin xmax ymin ymax]); pause(0.001);
end % function PlotXY





