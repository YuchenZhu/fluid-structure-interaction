clc; 
clear all;
close all;

%% Input data
N  = 20;        % number of cells

dt   = 0.1;     % time step
tend = 10;      % total simulation time

fmt = 'b-';     % linecolor for plot of final solution

show_sol = 1;   % (=1) show intermediate solutions during simulation

method = 1;     % method = 1 : exact at t(n+1)
                % method = 2 : exact at t(n+1/2)
                % method = 3 : DGCL
if (method == 1)
  disp('Mesh velocities are EXACT at t(n+1)');
elseif (method == 2)
  disp('Mesh velocities are EXACT at t(n+1/2)');
else
  disp('Mesh velocities satisfy DGCL');
end

%% Initialisation

dx = zeros(1,N);    % cell volumes
xi = zeros(1,N+1);  % face centers
x  = zeros(1,N);    % cell centers
u  = ones(1,N);     % solution



asdsda = (3:9)
xi  = (0:N)%/N;              % face centers
xi0 = (0:N)/N;              % face centers at t=0
x   = (1:N)/N - 0.5/N;      % cell centers
dx  = xi(2:N+1) - xi(1:N);  % cell volumes

asd=[0.5198
0.4757
0.2392
0.1191
0.0593];
q1=0.9152/0.4757;q1=q1^(-1);
q2=0.4757/0.2392;q2=q2^(-1);
q3=0.2392/0.1191;q3=q3^(-1);
q4=0.1191/0.0593;q4=q4^(-1);


x1=0.00022954495581;
x2=0.00170807576105;
x3=x1/x2;
0.134388041235883

xd=[1 2 3]
yd=[0.6 0.2 0.2]
loglog(xd,yd,'blacko-','MarkerFaceColor',[1,0,0]);grid on; axis([10^0 10^2 10^-5 10^0]); 