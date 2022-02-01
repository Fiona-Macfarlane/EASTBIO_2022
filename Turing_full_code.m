% This code was written and developed by F R Macfarlane 4/02/2021
%% Code to solve system of two PDES for two interacting chemicals u and v
% The system exhbits Turing patterns under parameter setting of code.
% 2D simulation 

% PDEs:
%   u_t = Du (u_xx + u_yy) +f(u,v)
%   v_t = Dv (v_xx + v_yy) +g(u,v)
% where
%   u_t, v_t = derivative of u, v with respect to time
%   u_xx, v_xx = second derivative of u, v with respect to x
%   u_yy, v_yy = second derivative of u, v with respect to y
% and the interaction functions are
%   f(u,v)=alphau-u +u^2 v
%   g(u,v)=alphav -u^2 v

%% Section 0: Clear all presets and set defaults (optional)
clear                       % Clear all current data within MATLAB
close all                   % Close all files open within MATLAB


set(0,'DefaultTextFontSize', 20);   % optional: Set default font size on figures
set(0, 'defaultFigureUnits', 'normalized');  % optional: set the positional units of a figure to be normalised (1=length/width of screen)
% next two lines set the deafult size of figures either medium or full
% screen ( the original default setting is quite small)
set(0, 'defaultFigurePosition', [0.2    0.2    0.7    0.7]) % medium
% set(0, 'defaultFigurePosition', [0    0    1    1]) % fullscreen
%% Section 1: Define domain
Nx=201;                     % Number of grid spaces in x direction
Ny=Nx;                      % Number of grid spaces in y direction (must be equal for the code to rrun as set up)
xmin=0;                     % Minimum x value
xmax=1;                     % Maximum x value
x=linspace(xmin,xmax,Nx);   % Set x grid
ymin=0;                     % Minimum y value
ymax=1;                     % Maximum y value
y=linspace(ymin,ymax,Ny);   % Set y grid
dx=(xmax-xmin)/(Nx-1);      % Find grid spacing in x direction
dy=(ymax-ymin)/(Ny-1);      % Find grid spacing in y direction
Grid=meshgrid(x,y);         % Define grid in x and y directions

dt=0.001;                   % Set time-step of simulations (note has to be >=0.001 for numerics to be stable)
T=200;                      % Set final time-step of simulations 

%% Section 2: Set initial conditions
u0=1;                       % Set the initial density (steady state value) for u
v0=0.9;                     % Set the initial density (steady state value) for v
ICScal=0.001;               % Set the scale of the perturbation
u=u0-ICScal+(2*ICScal)*rand(Nx,Ny); % Set u IC to be ± random perturbation of SS
v=v0-ICScal+(2*ICScal)*rand(Nx,Ny); % Set v IC to be ± random perturbation of SS

%% Optional: Set up vectors to store data over time
% This allows us to see the values of u and v over a range of time-points
% once the simulation has finished rather than just the final values

ustore={};                  % Set an empty MATLAB cell to store data for u over time
vstore={};                  % Set an empty MATLAB cell to store data for u over time
ustore{1}=u;                % Set initial data of u to be first stored value
vstore{1}=v;                % Set initial data of v to be first stored value
inds=1;                     % Set up an index to keep track of how many time-steps of data are stored

%% Section 3: Parameter settings
Du = 0.1*dt;    % Diffusion rate of chemical u (scaled by time-step length)
Dv = 4*dt;      % Diffusion rate of chemical v (scaled by time-step length)
alphau=0.1;     % Production rate of chemical u
alphav=0.9;     % Production rate of chemical v

%% Section 4: Run simulation
% Set a for loop for each time-step, where time=t*dt and final time is T
% That is the code within the loop will update at each time-step
for t=1:(T/dt)
    
    uold=u; % Set the concentration of u to be uold (the value of u at previous time-step)
    vold=v; % Set the concentration of v to be vold (the value of u at previous time-step)
    
    %% Section 4.1: Update interaction functions f and g
    % Set new variable ui, vi as the values of u and v on the domain excluding the
    % boundaries x=1,Nx, y=1,Ny.
    ui=uold(2:Nx-1,2:Ny-1);
    vi=vold(2:Nx-1,2:Ny-1);
    % Set the interaction functions, recall:
    %   f(u,v)=alphau-u +u^2 v
    %   g(u,v)=alphav -u^2 v
    K=1; % A potential scaling factor of our interaction terms
    f=K.*(alphau-ui+(ui.*ui.*vi));
    g=K.*(alphav-(ui.*ui.*vi));
    
    %% Section 4.2: Update diffusion terms
    % Use the discretised second derivative (finite difference)
    % approximation in space
    %   n(t,x,y)_xx \approx (n(t,x+1,y)-2 n(t,x,y) +n(t,x-1,y))/dx^2
    %   n(t,x,y)_yy \approx (n(t,x,y+1)-2 n(t,x,y) +n(t,x,y-1))/dy^2
    % We also note that the approximation in time is
    %   n(t,x,y)_t \approx (n(t+1,x,y)-n(t,x,y))/dt
    % So we multiply everything by dt so that we can evaluate n(t+1,x,y)

    % Define diffusion terms in x direction at all positions except those on the boundaries
    LapXu=Du.*(dt/(dx*dx)).*(uold(1:Nx-2,2:Ny-1)+uold(3:Nx,2:Ny-1)-2*uold(2:Nx-1,2:Ny-1));
    LapXv=Dv.*(dt/(dx*dx)).*(vold(1:Nx-2,2:Ny-1)+vold(3:Nx,2:Ny-1)-2*vold(2:Nx-1,2:Ny-1));
    
    % Define diffusion terms in y direction at all positions except those on the boundaries
    LapYu=Du.*(dt/(dy*dy)).*(uold(2:Nx-1,1:Ny-2)+uold(2:Nx-1,3:Ny)-2*uold(2:Nx-1,2:Ny-1));
    LapYv=Dv.*(dt/(dy*dy)).*(vold(2:Nx-1,1:Ny-2)+vold(2:Nx-1,3:Ny)-2*vold(2:Nx-1,2:Ny-1));
    
    %% Section 4.3: Update u and v using finte difference approximation
    % new values = old values + any change from interaction functions
    % (times dt because of approximation of time deriavtive) + any change
    % from diffusion
    
    % Update vectors u and v at all positions except those on the boundaries
    u(2:Nx-1,2:Ny-1)=uold(2:Nx-1,2:Ny-1)+dt.*f+LapXu+LapYu;
    v(2:Nx-1,2:Ny-1)=vold(2:Nx-1,2:Ny-1)+dt.*g+LapXv+LapYv;
    
    %% Section 4.4: Boundary conditions
    % Impose Zero Flux boundary conditions by setting values of u and v at boundaries
    % to be the same as those next to them  
    u(1,2:Ny-1)=u(2,2:Ny-1);
    u(Nx,2:Ny-1)=u(Nx-1,2:Ny-1);
    v(1,2:Ny-1)=v(2,2:Ny-1);
    v(Nx,2:Ny-1)=v(Nx-1,2:Ny-1);
    
    u(2:Nx-1,1)=u(2:Nx-1,2);
    u(2:Nx-1,Ny)=u(2:Nx-1,Ny-1);
    v(2:Nx-1,1)=v(2:Nx-1,2);
    v(2:Nx-1,Ny)=v(2:Nx-1,Ny-1);
    
    u(1,1)=u(2,2);
    u(Nx,Ny)=u(Nx-1,Ny-1);
    u(Nx,1)=u(Nx-1,2);
    u(1,Ny)=u(2,Ny-1);
    
    v(1,1)=v(2,2);
    v(Nx,Ny)=v(Nx-1,Ny-1);
    v(Nx,1)=v(Nx-1,2);
    v(1,Ny)=v(2,Ny-1);
    
    %% Optional: plot results and store data every 1/dt time-steps
    % Note, you could plot every time-step but this would slow down runs,
    % and changes would be minimal
    
    if mod(t,1/dt)==0 % if our time is a multiple of 1/dt run the plotting code
        
        % Plot concentration of u
        subplot(1,2,1)                      % We want two plots side by side subplot(number of rows of plots, number of columns of plots, which plot this is)
        surf(x,y,u,'EdgeColor','none')      % Surf gives us a surface plot of u on the domain x,y    
        axis square                         % Optional: this ensures a and y axis are same length in displayed plot   
        colormap jet;                       % This determines the colors used on the surf
        view([0,90]);                       % This decides the orientation of our plot - this option gives a top down view of results
        xlabel('x');                        % Label for x axis
        ylabel('y');                        % Label for y axis
        title(sprintf('Concentration of u at t = %d ', t*dt)); % Title of plot that updates for every time-step t
        hcb=colorbar;                       % Adds a colorbar for the concentration levels
        title(hcb,'Concentration of u');    % Adds label to colorbar
        box on                              % Optional: box around axis (plotting preference)
        
        % Plot concentration of v
        subplot(1,2,2)                      % We want two plots side by side subplot(number of rows of plots, number of columns of plots, which plot this is)
        surf(x,y,v,'EdgeColor','none')      % Surf gives us a surface plot of u on the domain x,y    
        axis square                         % Optional: this ensures a and y axis are same length in displayed plot   
        colormap jet;                       % This determines the colors used on the surf
        view([0,90]);                       % This decides the orientation of our plot - this option gives a top down view of results
        xlabel('x');                        % Label for x axis
        ylabel('y');                        % Label for y axis
        title(sprintf('Concentration of v at t = %d ', t*dt)); % Title of plot that updates for every time-step t
        hcb=colorbar;                       % Adds a colorbar for the concentration levels
        title(hcb,'Concentration of v');    % Adds label to colorbar
        box on                              % Optional: box around axis (plotting preference)
        
        
        drawnow                             % This ensures the figure will update throughout the simulation and not just at the end
        
        inds=inds+1;                        % Add one to the storage index
        ustore{inds}=u;                     % Add values of u to the storage vector
        vstore{inds}=v;                     % Add values of v to the storage vector
    end % end plotting if loop
end % end of time loop

%% Optional: Save data and final figure to the filenames in purple
save('Turing_Data.mat')
savefig('Turing_Figure.fig')

% END OF FILE %
