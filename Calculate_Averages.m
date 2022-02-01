%% File to find average and standard deviation between runs of the agent-based model of disease spread
% % This code was written by Dr Fiona R Macfarlane (Jan 2022)

% The file ABM_RunSingle must have been run more than once with the same parameter setting with data from
% each run saved in different files named Data*.mat *=1,..,n (Section 7) 
                                 
%% Section 0: Defaults (same as before)

clc                                         % clear previous commands
clear all                                   % clear previously stored information
close all                                   % close all open figures

set(0,'defaultFigureUnits','normalized');   % set default figure sizes to be normalised values
set(0,'defaultFigurePosition',[0 0 1 1]);   % set default figure size to be full screen
set(groot,'defaulttextinterpreter','tex');  % set interpreter of figures to allow for mathematical symbols  
set(0,'DefaultTextFontSize',25);            % set default font size for labels of figures
set(0,'DefaultAxesFontSize',20);            % set default font size for axes labels

%% Section 1: Load data from each file

NumFiles=10;                                % set the number of runs you want to compare

for a = 1:NumFiles                          % for each file you wish to compare
    load(sprintf('Data%d.mat',a),'NS','NI','NR','Final_Time') % load the required data from each file
    SVal{a} = NS;                           % save the number of susceptible agents over time for each file
    IVal{a} = NI;                           % save the number of infected/infective agents over time for each file
    RVal{a}= NR;                            % save the number of recovered agents over time for each file
    clear NS NI NR                          % clear the values to ensure they do not clash with next file
end

for b=1:Final_Time                          % for each time-step recorded
    for c=1:NumFiles                        % for each run you want to average
    Sus(c)=SVal{c}(b);                      % find the number of susceptible agents at time b in data set c
    Inf(c)=IVal{c}(b);                      % find the number of infected/infective agents at time b in data set c
    Rec(c)=RVal{c}(b);                      % find the number of recovered agents at time b in data set c
    end
    SMean(b)=mean(Sus);                     % find mean number of susceptible agents at time b across all runs of simulation
    SStd(b)=std(Sus);                       % find standard deviation of susceptible agents at time b across all runs of simulation
    IMean(b)=mean(Inf);                     % find mean number of infected/infective agents at time b across all runs of simulation
    IStd(b)=std(Inf);                       % find standard deviation of infected/infective agents at time b across all runs of simulation
    RMean(b)=mean(Rec);                     % find mean number of recovered agents at time b across all runs of simulation
    RStd(b)=std(Rec);                       % find standard deviation of recovered agents at time b across all runs of simulation    
end

%% Section 2: Plot the averages and standard deviaition

t=linspace(1,Final_Time,Final_Time);                        % set a time vector (for plotting standard deviation)
hold on  
plot(SMean,'b','LineWidth',3)                               % plot the mean number of susceptible agents over time as a blue line (b)
plot(IMean,'r','LineWidth',3)                               % plot the mean number of infected/infective agents over time as a red line (r)
plot(RMean,'k','LineWidth',3)                               % plot the mean number of recovered agents over time as a black line (k)
fill([t, fliplr(t)], [SMean + SStd, fliplr(SMean - SStd)], 'b','facealpha',.2);     % plot the standard deviation for susceptible agents as a blue transparent area
fill([t, fliplr(t)], [IMean + IStd, fliplr(IMean - IStd)], 'r','facealpha',.2);     % plot the standard deviation for infected/infective agents as a red transparent area
fill([t, fliplr(t)], [RMean + RStd, fliplr(RMean - RStd)], 'k','facealpha',.2);     % plot the standard deviation for recovered agents as a black transparent area
hold off  
% set legend for figure
legend('Mean number of susceptible agents','Mean number of infected agents','Mean number of recovered agents','Mean \pm SD susceptible agents','Mean \pm SD infected agents','Mean \pm SD recovered agents')
axis([1 Final_Time 0 600])                                  % set axis of plot to ensure all values can be seen
title(['Average number in each population over time with standard deviations, averaged over ',num2str(NumFiles),' runs'])    % set title of subplot
xlabel('Time')                                              % set x axis label
ylabel('Number in population')                              % set y axis label
box on                                                      % ensure box around figure

%saveas(gcf,'Figure_Averages.png')                      	 % optional to save the figure as a png file (to open outside matlab)
%savefig('Figure_Averages.fig')                              % optinoal to save the figure as a matlab figure that can then be opened/edited in matlab

% End of File %