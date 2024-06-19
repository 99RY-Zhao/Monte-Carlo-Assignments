%% EXERCISE 1 - MARKOV CHAIN MONTE CARLO FOR COAL MINE ACCIDENTS
close all;
clearvars;
clc;

format compact;
load coal_mine_disasters.mat;
addpath(genpath("functions"));

% ESTIMATION PROCESS HYPERPARAMETER DEFINITION.
limits = [1658 1980]; % Define time interval.
d = 6;                % Define break count.
psi = 8;              % Define theta distrib.
rho = 0.05;           % Define t marginal dist.
r = "individual";     % MH proposal kernel ("global" or "individual").

% SIMULATION LENGTH AND BURN-IN
burnin = 1e4;        % Burn-in
obsrvs = 1e5;        % Observations
N = burnin + obsrvs; % Total time

% Perform the MCMC algorithm
[history,~] = coalmine(T,d,psi,rho,r,N,true);

% Extrapolate the history of all elements
breakhist = vertcat(history{:,3});
lambdahist = vertcat(history{:,2});
thetahist = vertcat(history{:,1});

% Estimate the quantities
t = mean(breakhist(burnin+1:end,:),1);
lambda = mean(lambdahist(burnin+1:end,:),1);
theta = mean(thetahist(burnin+1:end,:),1);



%% OBSERVING VARIANCE OF ALL ELEMENTS
fprintf("Var(lambda) = ");
fprintf("%e, ",var(lambdahist(burnin+1:end,:),1));

fprintf("\nVar(t) = ");
fprintf("%e, ",var(breakhist(burnin+1:end,:),1));

fprintf("\nVar(theta) = %e\n",var(thetahist(burnin+1:end,:),1));



%% PLOTTING THE HISTORY OF BREAKPOINTS
figure(1), hold on;
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
title("Evolution of breakpoints");

% Set axes limits and labels
ylim(limits);
ylabel("Year");

xlim([1,N]);
xlabel("Epoch");



%% PLOTTING THE HISTORY OF LAMBDA
figure(2), hold on;

for i=1:d
    plot(1:N,lambdahist(:,i));
end

% Calculate the y-scale of the plot
endpoints = [min(min(lambdahist)) max(max(lambdahist))];
plotwidth = endpoints(2) - endpoints(1);
endpoints = [endpoints(1)-0.1*plotwidth endpoints(2)+0.1*plotwidth];

% Display end of burn-in
plot([burnin burnin],endpoints,"Color","red");

% Adjust properties
title("Evolution of intensities");

% Set axes limits and labels
ylim(endpoints);
ylabel("lambda");

xlim([1,N]);
xlabel("Epoch");



%% PLOTTING THE HISTORY OF THETA
figure(3), hold on;
plot(1:N,thetahist);

% Calculate the y-scale of the plot
endpoints = [min(thetahist) max(thetahist)];
plotwidth = endpoints(2) - endpoints(1);
endpoints = [endpoints(1)-0.1*plotwidth endpoints(2)+0.1*plotwidth];

% Display end of burn-in
plot([burnin burnin],endpoints,"Color","red");

% Adjust properties
title("Evolution of theta");

% Set axes limits and labels
ylim(endpoints);
ylabel("theta");

xlim([1,N]);
xlabel("Epoch");



%% PLOTTING THE FITTED INTENSITIES TO THE CURVE
figure(4), hold on;

% Plot the curve form the data
plot(T,1:length(T));

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
title("Cumulative number of accidents");

% Set axes limits and labels
ylim([0 length(T)]);
ylabel("Accidents");

xlim([T(1) T(end)]);
xlabel("Year");



%% OBSERVE EFFECT ON POSTERIOR VARIANCE OF PSI AND RHO CHOICE
d = 6;            % Define break count.
r = "individual"; % MH proposal kernel ("global" or "individual").

% SIMULATION LENGTH AND BURN-IN
burnin = 1e4;        % Burn-in
obsrvs = 1e4;        % Observations
N = burnin + obsrvs; % Total time


psigrid = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40];
rhogrid = [0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];

%% Perform the tests
results = cell(length(psigrid),length(rhogrid),3);
for psi=1:length(psigrid)
    for rho=1:length(rhogrid)
        [history,~] = coalmine(T,d,psigrid(psi),rhogrid(rho),r,N);

        results{psi,rho,3} = var(vertcat(history{:,3}),1);
        results{psi,rho,2} = var(vertcat(history{:,2}),1);
        results{psi,rho,1} = var(vertcat(history{:,1}),1);
    end
end

%% Map and average the results to a matrix
thetamap = cell2mat(results(:,:,1));
lambdamap = cellfun(@mean,results(:,:,2));
breakmap = cellfun(@mean,results(:,:,3));

% Plot everything
figure(5);
surf(rhogrid,psigrid,thetamap);
title("Variance of theta with respect to rho and psi");
xlabel("rho");
ylabel("psi");

figure(6);
surf(rhogrid,psigrid,lambdamap);
title("Average variance of lambda with respect to rho and psi");
xlabel("rho");
ylabel("psi");

figure(7);
surf(rhogrid,psigrid,breakmap);
title("Average variance of t with respect to rho and psi");
xlabel("rho");
ylabel("psi");



%% EXERCISE 2 - PARAMETRIC BOOTSTRAP FOR THE 100-YEAR ATLANTIC WAVE
close all;
clearvars;
clc;

load("atlantic.txt");

% Define inverse function
F_inv = @(u,mu,beta) mu - beta*log(-log(u));

n = length(atlantic);
B = 1000;

[beta_hat, mu_hat] = est_gumbel(atlantic);
boot = zeros(2,B);
for b = 1:B
    y_boot = F_inv(rand(n,1),mu_hat,beta_hat);
    [beta_hat_, mu_hat_] = est_gumbel(y_boot);
    boot(1,b) = beta_hat_;
    boot(2,b) = mu_hat_;
end

delta_beta = sort(boot(1,:) - beta_hat);
delta_mu = sort(boot(2,:) - mu_hat);
alpha = 0.05;

LB_beta = beta_hat - delta_beta(ceil((1 - alpha/2)*B));
UB_beta = beta_hat - delta_beta(ceil(alpha*B/2));

LB_mu = mu_hat - delta_mu(ceil((1 - alpha/2)*B));
UB_mu = mu_hat - delta_mu(ceil(B*alpha/2));

beta_LB_UB = [LB_beta UB_beta (UB_beta - LB_beta)];
mu_LB_UB = [LB_mu UB_mu (UB_mu - LB_mu)];

T = 3 * 14 * 100;
wave_height = zeros(1,B);
for b = 1:B
    wave_height(1,b) = F_inv(1 - 1/T,boot(2,b),boot(1,b));
end

wave_hat = F_inv(1 - 1/T,mu_hat,beta_hat);
delta_wave = sort(wave_height - wave_hat);

UB_wave = wave_hat - delta_wave(ceil(alpha*B));




beta_LB_UB
mu_LB_UB
UB_wave

% PLOTS !!