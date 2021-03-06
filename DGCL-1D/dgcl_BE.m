%% Input data
N  = 20;        % number of cells







method = 1;     % method = 1 : exact at t(n+1)
                % method = 2 : exact at t(n+1/2)
                % method = 3 : DGCL
k = 0;
timesteps = 0.1*2.^( -1* [ 0 1 2 3 4 5 6 ] ) ;
for dt = timesteps
    tend = 10;      % total simulation time
    if dt == 0.1
        fmt = 'b-';     % linecolor for plot of final solution
    elseif dt ==0.05
        fmt = '-r';
    elseif dt == 0.025
        fmt = '-g';
    end
    
    show_sol = 0;   % (=1) show intermediate solutions during simulation
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

xi  = (0:N)/N;              % face centers
xi0 = (0:N)/N;              % face centers at t=0
x   = (1:N)/N - 0.5/N;      % cell centers
dx  = xi(2:N+1) - xi(1:N);  % cell volumes

% Plot initial solution
figure(1);
hold off;
axis([0 1 0.98 1.02]);
plot(x,u,'-x');
hold off;

L = zeros(N);   % Discretization matrix

%% Simulation loop
j=0;
for t=dt:dt:tend
  dx_tn = dx;  % store cell volume at tn
  xi_tn = xi;  % store face center at tn
  u_tn  = u;   % store solution at tn
  
  xi     = xi0 + sin(2*pi*t) * sin(2*pi*xi0) / N;   % face centers at tn+1
  dx     = xi(2:N+1)-xi(1:N);                       % cell volumes at tn+1
  x      = xi(1:N) + dx/2;                          % cell centers at tn+1
  
  dxidt_exnp1   = 2*pi*cos(2*pi*(t))*sin(2*pi*xi0) / N;         % exact face velocity at tn+1
  dxidt_exnp1_2 = 2*pi*cos(2*pi*(t-dt/2))*sin(2*pi*xi0) / N;    % exact face velocity at tn+1/2
  dxidt_dgcl    = (xi - xi_tn) / dt;                            % face velocity satisfying D-GCL
  
  if (method == 1)
    dxidt = dxidt_exnp1;
  elseif (method ==2)
    dxidt = dxidt_exnp1_2;
  else
    dxidt = dxidt_dgcl;
  end
    
  % Setting up system to solve; for each cell:
  %
  %   (u*dx)^(n+1) - (u*dx)^n      
  %   ------------------------- - [ u_leftface * dxidt_leftface * (-1) + u_rightface * dxidt_rightface * (+1) ]^(n+1) = 0
  %             dt                 
  %
  % which can be written for the total system as:
  %
  % L u^(n+1) = (u*dx)^n
  %
  % We use a simple avarage to compute the solution at a cell face:
  %   u_face = 0.5 * (u_leftcell + u_rightcell)
  %
  % Internal cells
  for i=2:N-1
    L(i,i-1:i+1) = [0 dx(i) 0] - dt * 0.5*[-dxidt(i) dxidt(i+1)-dxidt(i) dxidt(i+1)];
  end
  % Boundary cells
  L(1,1:2)   = [dx(1) 0] - dt * 0.5*[dxidt(2)-2*dxidt(1) dxidt(2)];
  L(N,N-1:N) = [0 dx(N)] - dt * 0.5*[-dxidt(N) 2*dxidt(N+1)-dxidt(N)];

  % solve system
  u = (L\(dx_tn.*u_tn)')';

  % Show intermediate solution
  if (show_sol)
    figure(1);
    hold off;
    plot(x,u,'-x');
    axis([0 1 0.8 1.2]);
%    pause
  end
%Calculate first error
  e1 = norm(u - ones(1,N))/sqrt(N);
  j=j+1 ; E1(k+1, j) = e1 ;
end
%calculate the second error
k=k+1;
e2 = sqrt(sum(E1(k,:).^2)/length(0:dt:tend));%second error for a dt 
E2(k) = e2;
EE1(k) = e1;
%% Plot final solution for a time step
figure(2);
hold on;
plot(x,u,fmt);
title(['Solution at t=' num2str(tend)]);
ylabel('Solution');
xlabel('x');
Legend{k} = strcat('\Delta t=', num2str(dt));legend(Legend);
end
%% Plotting for error1 and 2
figure(3);
loglog(timesteps, EE1, '-rp',timesteps,E2,'-bp');
grid on;
xlabel('\Delta t'); ylabel('Error');
legend('error 1','error2')

%finding average slope for error 1 and 2
for i=1:6
    slope1(i) = log(EE1(i+1)/EE1(i))/log(timesteps(i+1)/timesteps(i));
     slope2(i) = log(E2(i+1)/E2(i))/log(timesteps(i+1)/timesteps(i));
end
avg_slope1 = mean(slope1);
avg_slope2 = mean(slope2);

