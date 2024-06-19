%% BEHAVIOR WRT BREAKPOINTS
close all;
clearvars;
clc;

format compact;
load coal_mine_disasters.mat;
addpath(genpath("functions"));

limits = [1658 1980]; % Define time interval.
psi = 8;              % Define theta distrib.
i = 0.05;           % Define t marginal dist.
r = "individual";     % MH proposal kernel ("global" or "individual").

% SIMULATION LENGTH AND BURN-IN
burnin = 1e4;        % Burn-in
obsrvs = 1e5;        % Observations
N = burnin + obsrvs; % Total time

% Perform the MCMC algorithm
for d=2:7
    history = coalmine(T,d,psi,i,r,N,false);

    % Extrapolate the history of all elements
    breakhist = vertcat(history{:,3});
    lambdahist = vertcat(history{:,2});
    thetahist = vertcat(history{:,1});

    % Estimate the quantities
    t = mean(breakhist(burnin+1:end,:),1);
    lambda = mean(lambdahist(burnin+1:end,:),1);
    theta = mean(thetahist(burnin+1:end,:),1);
 
    fprintf("%d    %e\n",d,var(thetahist(burnin+1:end,:),1));

    figure, hold on, subplot(1,2,1), hold on;
    breakhist = vertcat(history{:,3});

    % Start with the boundaries
    plot(1:N,breakhist(:,1),"Color","b");
    plot(1:N,breakhist(:,end),"Color","b");

    for i=2:(size(breakhist,2)-1)
        plot(1:N,breakhist(:,i));
    end

    % Display end of burn-in
    plot([burnin burnin],limits,"Color","red");

    % Adjust properties
    title("Evolution of breakpoints for d=" + d);

    % Set axes limits and labels
    ylim(limits);
    ylabel("Year");

    xlim([1,N]);
    xlabel("Epoch");

    subplot(1,2,2), hold on;
    
    % Plot the curve form the data
    plot(T,1:length(T),"r");
    
    % Fit the intensities
    previous = 0;
    
    for i=1:d
        period = t(i+1) - t(i);
        next = previous + lambda(i)*period;
    
        % Draw the endpoint
        plot(t(i+1),next,"b:","MarkerSize",25);
    
        % Draw the line
        plot([t(i) t(i+1)],[previous next],"b");
    
        previous = next;
    end
    
    % Adjust properties
    title("Cumulative number of accidents (d=" + d + ")");
    
    % Set axes limits and labels
    ylim([0 length(T)]);
    ylabel("Accidents");
    
    xlim([T(1) T(end)]);
    xlabel("Year");
end

%% ACCEPTANCE
close all;
clearvars;
clc;

format compact;
load coal_mine_disasters.mat;
addpath(genpath("functions"));

limits = [1658 1980]; % Define time interval.
d = 6;                % Breakpoint count.
psi = 8;              % Define theta distrib.
r = "individual";     % MH proposal kernel ("global" or "individual").

% SIMULATION LENGTH AND BURN-IN
burnin = 1e4;        % Burn-in
obsrvs = 1e4;        % Observations
N = burnin + obsrvs; % Total time

% Perform the MCMC algorithm
rhogrid = [0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
acceptances = zeros(size(rhogrid));
for i=1:length(rhogrid)
    [~,acc] = coalmine(T,d,psi,rhogrid(i),r,N,false);
    acceptances(i) = acc;
end

figure, hold on

for i=1:length(rhogrid)
    plot(rhogrid(i),acceptances(i),"b:","MarkerSize",25);

    if i > 1
        plot([rhogrid(i-1) rhogrid(i)],[acceptances(i-1) acceptances(i)],"b");
    end
end

title("Acceptance rate with respect to rho");

ylim([0 1]);
ylabel("Acceptance");

xlim([min(rhogrid) max(rhogrid)]);
xlabel("rho");