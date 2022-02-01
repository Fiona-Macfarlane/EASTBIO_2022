%% File to run agent-based model of disease spread
% % This code was written by Dr Fiona R Macfarlane (Jan 2022)

% This file runs a simulation of a population of people containing:
% Population S: Those susceptible to some virus
% Population I: Those infected by a virus and able to pass on the virus
% Population R: Those recovered from the virus that can no longer get it

%% Section 0: Defaults 

clc                                         % clear previous commands
clear all                                   % clear previously stored information
close all                                   % close all open figures

set(0,'defaultFigureUnits','normalized');   % set default figure sizes to be normalised values
set(0,'defaultFigurePosition',[0 0 1 1]);   % set default figure size to be full screen
set(groot,'defaulttextinterpreter','tex');  % set interpreter of figures to allow for mathematical symbols  
set(0,'DefaultTextFontSize',25);            % set default font size for labels of figures
set(0,'DefaultAxesFontSize',20);            % set default font size for axes labels
set(0,'DefaultUicontrolFontsize',18);       % set font size for pop up menu

%% Section 1: Set up the spatial domain

xm=0;                                       % minimum x value on grid
xM=100;                                     % maximum x value on grid
ym=0;                                       % minimum y value on grid
yM=100;                                     % maximum y value on grid

dx=1;                                       % length of each grid-space in x direction
dy=1;                                       % length of each grid-space in y direction

Nx=xM/dx;                                   % calculate number of positions in x direction
Ny=xM/dy;                                   % calculate number of positions in y direction

x=xm:dx:xM;                                 % define x grid 
y=ym:dy:yM;                                 % define y grid 


%% Section 3: Set up total number and initial positions of the population

N=400;                                      % set number of agents (people) in simulation

nx=zeros(1,N);                              % set up vector to store the x positions of each agent
ny=zeros(1,N);                              % set up vector to store the y positions of each agent

% Set up a menu so the user can choose whether or not to reset the initial
% condition or use a presaved value
info = menu('Do you want to reset the initial spatial positions of the population? Note, if this is the first time you are running file select yes', 'Yes', 'No');

if info == 1                                % if you are resetting the initial condition
    for t=1:N                               % start FOR loop for each agent:
        nx(t)=randi(Nx);                    % pick a random initial x position
        ny(t)=randi(Nx);                    % pick a random initial y position
        
        A=find(nx==nx(t));                  % find how many agents have this same x position
        B=find(ny==ny(t));                  % find how many agents have this same y position
        C=intersect(A,B);                   % find how mant agents have this same x AND y position
        
        while length(C)>1                   % start WHILE loop for cases where there is more than 1 agent on this position
            nx(t)=randi(Nx);                % pick a new x position
            ny(t)=randi(Nx);                % pick a new y position
            
            A=find(nx==nx(t));              % find how many agents have this same new x position
            B=find(ny==ny(t));              % find how many agents have this same new y position
            C=intersect(A,B);               % find how mant agents have this same new x AND y position
        end                                 % once you find an unoccupied position then the loop ends
    end                                     % once each agent has been placed without overlaps the FOR loop ends
    
    save('Initial_condition.mat','nx','ny');% save the initial condition so that we can use this next time if required

else                                        % if you are not resetting the intial condition
    
    load('Initial_condition.mat','nx','ny');% load the presaved initial condition

end                                         % end if condition (line 48)

%% Section 4: Set up the initial infection in the population

Infection=zeros(1,N);                       % set up a vector to store the infection status of each agent in the population (0 if uninfected, 1 if infected)                            

Time_Infected=zeros(1,N);                   % set up a vector to store the time since infection for each agent

Initial_No_Infected=10;                     % set the initial number of people who have the infection (minimum 1)

Vaccinated=0;                               % set the initial number of people who are recovered immediately (e.g. vaccinated)

Infection(1:Initial_No_Infected)=1;         % set the infection status of the desired number of people to be 1 (infected)

if Vaccinated>0
Infection(Initial_No_Infected+1:Initial_No_Infected+1+Vaccinated)=2;        % set the infection status of the desired number of vaccinated people to be 2 (recovered)
end

NS=length(find(Infection==0));              % calculate the initial number of susceptible agents
NI=length(find(Infection==1));              % calculate the initial number of infected agents
NR=length(find(Infection==2));              % calculate the initial number of recovered agents

nxstore{1}=nx;                              % set up 'cell' to store x positions of population over time
nystore{1}=ny;                              % set up 'cell' to store y positions of population over time
Infectionstore{1}=Infection;                % set up 'cell' to store infection status of population over time

%% Section 5: Set the parameter values of the simulation

Infection_Probability=1.0;                  % set the probability that a susceptible person will become infected by an infective person (value between 0 and 1)

Recovery_Time=50;                           % set the time in which an infected person is infective (real number greater than 1)

Final_Time=200;                             % set how many time-steps you want to run the simulations for (real number greater than or equal 1)

Radius=1;                                   % set the radius of infection (how close an infected person can infect people) (minimum 1, real values)

Mov=1;                                      % set probability of movement for each agent (value between 0 and 1)

%% Section 6: Run the time-loop simulations

for t=1:Final_Time                          % for each time-step t of the simulations run the following mechanisms
    
    %% Section 6a: Allow movement of population
    % Every agent has equal probability (Mov) of moving to the 4 neighbouring
    % positions at any time
    
    Md1=Mov/4;                              % define probability of moving left
    Md2=Mov/2;                              % define probability of moving right
    Md3=Mov*(3/4);                          % define probability of moving down
    Md4=Mov;                                % define probability of moving up
    
    Rad=rand(1,N);                          % generate one random number for each agent
    
    for a=1:N                               % for all people
        
        % The following lines check if space surrounding cell is free
        
        VX=find(nx==nx(a));                 % find how many agents have this x position
        VY=find(ny==ny(a));                 % find how many agents have this y position
        VXp=find(nx==nx(a)+1);              % find how many agents have this x position
        VYp=find(ny==ny(a)+1);              % find how many agents have this y position
        VXm=find(nx==nx(a)-1);              % find how many agents have this x position
        VYm=find(ny==ny(a)-1);              % find how many agents have this y position
        
        VU=intersect(VX,VYp);               % find how mant agents have this x AND y position
        VD=intersect(VX,VYm);               % find how mant agents have this x AND y position
        VL=intersect(VXm,VY);               % find how mant agents have this x AND y position
        VR=intersect(VXp,VY);               % find how mant agents have this x AND y position
        
        if Rad(a)<Md1 && length(VL)==0 && (nx(a)-1)>=1     % If random number is less than probability and the left position is free and the desired position is within boundary
        nx(a)=nx(a)-1;                      % Move left
        
        elseif Rad(a)>=Md1 && Rad(a)<Md2 && length(VR)==0  && (nx(a)+1)<=Nx  % If random number is less than probability and the right position is free and the desired position is within boundary
        nx(a)=nx(a)+1;                      % Move right
        
        elseif Rad(a)>=Md2 && Rad(a)<Md3 && length(VD)==0 && (ny(a)-1)>=1 % If random number is less than probability and the lower position is free and the desired position is within boundary
        ny(a)=ny(a)-1;                      % Move down
        
        elseif Rad(a)>=Md3 && Rad(a)<Md4 && length(VU)==0 && (ny(a)+1)<=Ny  % If random number is less than probability and the upper position is free and the desired position is within boundary
        ny(a)=ny(a)+1;                      % Move up
    end

    end
  
    
    %% Section 6b: Allow infection and recovery of population
    
    for l=1:N                               % for all agents in the population (agent l)
        for j=1:N                           % for all agents in the population (agent j to compare to agent l)
            
            A=rand;                         % set a random number for each agent
            
            % if the following are true then agent l becomes infected
                % the distance between the x positions of agent l and agent j is less than Radius*dx
                % the distance between the y positions of agent l and agent j is less than Radius*dy
                % agent l is uninfected
                % agent j is infected
                % the value of the random number A is less than the infection probability
            if round(abs(nx(l)-nx(j)))<=(Radius*dx) && round(abs(ny(l)-ny(j)))<=(Radius*dy) && Infection(l)==0 && Infection(j)==1 && Infection_Probability>=A
                
                Infection(l)=1;             % set the infection status of agent l to be 1 (infected)
                
            end                             % end loop (if loop line 161)
        end                                 % end loop comparing for agent j (for loop line 150)
        
        
        % If the agent is infected add 1 to the 'time infected' in vector Time_Infected
        if Infection(l)==1
            Time_Infected(l)=Time_Infected(l)+1;
        end
        
        % If the agent has been infected for longer than the set recovery time then the agent becomes recovered
        if Time_Infected(l)>Recovery_Time
            Infection(l)=2;
        end
    end
    
    %% Section 6c: Find spatial positions and total number of each population (for plotting)
    
    Sx=nx(find(Infection==0));              % find the x positions of suceptible population
    Sy=ny(find(Infection==0));              % find the y positions of suceptible population
    
    Ix=nx(find(Infection==1));              % find the x positions of infected/infective population
    Iy=ny(find(Infection==1));              % find the y positions of infected/infective population
    
    Rx=nx(find(Infection==2));              % find the x positions of recovered population
    Ry=ny(find(Infection==2));              % find the y positions of recovered population
    
    NS=[NS length(Sx)];                     % add the number of susceptible agents to the vector storing this over time
    NI=[NI length(Ix)];                     % add the number of infected/infective agents to the vector storing this over time
    NR=[NR length(Rx)];                     % add the number of recovered agents to the vector storing this over time
    
    %% Section 6d: Plot results each time-step and save data
    
    % figure:
    
    clf                                     % clear any previous data on figure (at each time-step)
    
    % subplot 1 displays the spatial positions of agents
    
    subplot(3,2,[1,2,3,4])                  % set up first subplot (define position of plot)
    hold on                                 % hold on allows you to plot multiple data on the same plot
    plot(Sx,Sy,'b.','MarkerSize',15)        % plot the spatial positions of susceptible agents as blue dots
    plot(Ix,Iy,'r.','MarkerSize',15)        % plot the spatial positions of infected/infective agents as red dots
    plot(Rx,Ry,'k.','MarkerSize',15)        % plot the spatial positions of recovered agents as blue dots
    hold off                                % end hold
    axis square                             % to enforce the plot to be square
    xticks([])                              % remove the x axis labels
    yticks([])                              % remove the y axis labels
    box on                                  % ensure the plot has box around it
    title(['Spatial positions at time=',num2str(t)])        % set title of plot that changes each time-step
    % set a string for key on the right hand side 
    str = {['\color{blue} \bullet Susceptible'],[],['\color{red} \bullet Infected/Infective'],[],['\color{black} \bullet Recovered']};
    % plot the key of the spatial positions plot in desired position
    annotation('textbox',[0.7 0.65 0.1 0.1],'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
    
    % subplot 2 will display the running total number in each population over time
    
    subplot(3,2,5.5)                        % set up subplot (define position of plot)
    hold on                                 % hold on function allows multiple data on same plot
    plot(NS,'b','LineWidth',3)              % plot the total number of susceptible agents over time as a blue line (b)
    plot(NI,'r','LineWidth',3)              % plot the total number of infected/infective agents over time as a red line (r)
    plot(NR,'k','LineWidth',3)              % plot the total number of susceptible agents over time as a black line (k)
    hold off                                % end hold function
    axis([1 Final_Time 0 600])              % set axis of plot to ensure all values can be seen
    % Set string for key for this subplot
    str = {['\color{blue} - Susceptible'],['\color{red} - Infected/Infective'],['\color{black} - Recovered']};
    % plot the key of the total numbers plot in desired position
    annotation('textbox',[0.7 0.2 0.1 0.1],'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
    title(['Total number in each population over time'])    % set title of subplot
    xlabel('Time')                          % set x axis label
    ylabel('Number in population')          % set y axis label
    box on                                  % ensure that plot has a box around it
    
    drawnow                                 % drawnow ensures that the plot is displayed after each time-step
    
    % Store data in the pre-assigned cells to allow for storage of data for future comparison 
    
    nxstore{t+1}=nx;                        % add the x positions of all agents for this time-step
    nystore{t+1}=ny;                        % add the y positions of all agents for this time-step
    Infectionstore{t+1}=Infection;          % add the infection status of all agents for this time-step
    
    
end                                         % end time-loop (if time < Final_Time then the loop will run again)

%% Section 7: Optional save data and figure

                        
save('Data1.mat');                                           % option 1: save all matlab data to file called Data (comment out if not using)
% save('Data1.mat','nxstore','nystore','Infectionstore','NS','NI','NR','Final_Time');    % option 2: save only the spatial positions, infection status and populations over time to matlab data to file called Data (uncomment if using)

%saveas(gcf,'Figure1.png')                                    % optional to save the figure as a png file (to open outside matlab)
%savefig('Figure1.fig')                                       % optinoal to save the figure as a matlab figure that can then be opened/edited in matlab

%% Section 8: Optional alert that file is finished
beep;pause(2);beep;                                         % the computer will beep twice in succession to alert user the code has finished                                   
% close all                                                 % optionally close everything once finished (uncomment if desired)

% End of File %
