%% PROBLEMS 3-6: ESTIMATION OF Cn WITH SIS(R)
close all;
clearvars;
clc;

addpath("functions");

N = 10000; % Number of particles for experiment
R = 10;    % Number of repetitions for each experiment
M = 24;    % Number of steps to simulate to.

% Cache results
if exist("results2.mat","file") == 2
    load("results2.mat");
else
    results = cell(1,M);
    for steps=1:M
        outcomes = zeros(R,3);
        for repetition=1:R
            [Xn,Yn,Wn] = sawsis2(N,steps,"naive",false);
            [Xs,Ys,Ws] = sawsis2(N,steps,"smart",false);
            [Xr,Yr,Cr] = sawsis2(N,steps,"smart",true);

            Cn = mean(Wn(end,:));
            Cs = mean(Ws(end,:));

            disp("    " + Cn + "/" + Cs + "/" + Cr);
            outcomes(repetition,:) = [Cn Cs Cr];
        end
    
        means = mean(outcomes,1);
        disp("[" + steps + "]: Average on " + R + " is " + means(1) + "/" + means(2) + "/" + means(3));
    
        results{steps} = outcomes;
    end

    save("results2.mat","results");
end

% Gather data on the experiment
data = -1 * ones(M,9);
for steps=1:M
    for repetition=1:R
        r = results{steps}(repetition,:);
        for value=1:3
            if data(steps,1+((value-1)*3)) == -1 || data(steps,1+((value-1)*3)) > r(value)
                data(steps,1+((value-1)*3)) = r(value);
            end
            if data(steps,value*3) == -1 || data(steps,value*3) < r(value)
                data(steps,value*3) = r(value);
            end
            data(steps,2+((value-1)*3)) = data(steps,2+((value-1)*3)) + r(value);
        end
    end
    data(steps,[2,5,8]) = data(steps,[2,5,8])/R;
end

% Make a nice plot
figure(1), hold on, legend show;
xlabel("Number of steps");
ylabel("log10(c_n)");

title("Estimation of c_n using the three discussed methods")
plot(1:24,log10(data(:,1)),":b","HandleVisibility","off")
plot(1:24,log10(data(:,2)),"b","DisplayName","SIS + trivial g")
plot(1:24,log10(data(:,3)),":b","HandleVisibility","off")
plot(1:24,log10(data(:,4)),":r","HandleVisibility","off")
plot(1:24,log10(data(:,5)),"r","DisplayName","SIS + smart g")
plot(1:24,log10(data(:,6)),":r","HandleVisibility","off")
plot(1:24,log10(data(:,7)),":g","HandleVisibility","off")
plot(1:24,log10(data(:,8)),"g","DisplayName","SISR")
plot(1:24,log10(data(:,9)),":g","HandleVisibility","off")

% Setup the linear regression
X = [ones(M-1,1) (1:M-1)' log(1:M-1)'];

% For each repetition of the experiment we get three parameters:
% IN THIS ORDER: [log(A_d) log(mu_d) (gamma_d-1)]
ThetaSIS  = zeros(R,3);
ThetaSISR = zeros(R,3);

for index=1:R
    TargetSIS  = zeros(M,1);
    TargetSISR = zeros(M,1);

    for steps=1:M
        TargetSIS(steps)  = log(results{steps}(index,2));
        TargetSISR(steps) = log(results{steps}(index,3));
    end

    ThetaSIS(index,:)  = transpose(X\TargetSIS(2:end));
    ThetaSISR(index,:) = transpose(X\TargetSISR(2:end));
end

% Adjust the estimates using their formula
ParamsSIS  = [exp(ThetaSIS(:,1))  exp(ThetaSIS(:,2))   1+ThetaSIS(:,3)];
ParamsSISR = [exp(ThetaSISR(:,1)) exp(ThetaSISR(:,2))  1+ThetaSISR(:,3)];

% Calculate variances
VarianceSIS  = var(ParamsSIS,1,1);
VarianceSISR = var(ParamsSISR,1,1);

% Now derive the actual measures
PSIS  = mean(ParamsSIS,1);
PSISR = mean(ParamsSISR,1);

display("A_2 = " + PSIS(1) + "/" + PSISR(1) + " Var: " + VarianceSIS(1) + "/" + VarianceSISR(1));
display("mu_2 = " + PSIS(2) + "/" + PSISR(2) + " Var: " + VarianceSIS(2) + "/" + VarianceSISR(2));
display("gamma_2 = " + PSIS(3) + "/" + PSISR(3) + " Var: " + VarianceSIS(3) + "/" + VarianceSISR(3));



%% PROBLEM 9: ESTIMATION OF PARAMETERS FOR 3 DIMENSIONS
close all;
clearvars;
clc;

addpath("functions");

N = 10000; % Number of particles for experiment
R = 10;    % Number of repetitions for each experiment
M = 24;    % Number of steps to simulate to.

% Cache results
if exist("results3.mat","file") == 2
    load("results3.mat");
else
    results = zeros(R,M);
    for steps=1:M
        outcomes = zeros(R,1);
        for repetition=1:R
            [X,Y,Z,C] = sawsis3(N,steps,"smart",true);
    
            disp("  |-> " + C);
            results(repetition,steps) = C;
        end

        disp("[" + steps + "]: Average is " + mean(results(:,steps)));    
    end

    save("results3.mat","results");
end

% Setup the linear regression
X = [ones(M-1,1) (1:M-1)' log(1:M-1)'];

% For each repetition of the experiment we get three parameters:
% IN THIS ORDER: [log(A_d) log(mu_d) (gamma_d-1)]
Theta  = zeros(R,3);

for index=1:R
    Theta(index,:)  = transpose(X\log(results(index,2:M)'));
end

Readings = [exp(Theta(:,1)) exp(Theta(:,2)) 1+Theta(:,3)];

ParVar      = var(Readings,1,1);
Parameters  = mean(Readings,1);

display("A_3 = " + Parameters(1) + " Var: " + ParVar(1));
display("mu_3 = " + Parameters(2) + " Var: " + ParVar(2));
display("gamma_3 = " + Parameters(3) + " Var: " + ParVar(3));

%% PROBLEM 10: SISR Implementation of Linear/Gaussian HMM.
close all;
clearvars;
clc;

load('population_2023.mat');

N = 1000;
n = 50;
A = 0.9; B = 3.9; C = 0.6; D = 0.99; G = 0.7; H = 1.2;

tau = zeros(1,n+1);
LB = zeros(1,n+1);
UB = zeros(1,n+1);

p = @(x, y) unifpdf(y, G*x, H*x);

part = unifrnd(C,D,N,1);  % initialization
w = p(part,Y(1)); % weighting
tau(1) = sum(part.*w)/sum(w); %estimation

[x_sort, i] = sort(part);
cum_ = cumsum(w(i))/sum(w);

I_L = find(cum_ >= 0.025, 1); % Index of lower 2.5% quantile
I_U = find(cum_ >= 0.975, 1); % Index of upper 2.5% quantile
LB(1) = x_sort(I_L);
UB(1) = x_sort(I_U);

ind = randsample(N,N,true,w);
part = part(ind);

for k = 2:n+1
    part = unifrnd(A,B,N,1).*part.*(1-part); % Mutation
    w = p(part, Y(k));
    tau(k) = sum(part.*w) / sum(w);

    [x_sort, i] = sort(part);
    cum_ = cumsum(w(i))/sum(w);

    I_L = find(cum_ >= 0.025, 1); % Index of lower 2.5% quantile
    I_U = find(cum_ >= 0.975, 1); % Index of upper 2.5% quantile
    LB(k) = x_sort(I_L);
    UB(k) = x_sort(I_U);

    ind =  randsample(N,N,true,w);
    part = part(ind);
end

oo = [LB;X';UB];

figure()
hold on
plot(1:n+1, tau,  '--d', 'LineWidth', 2,'Color','#77DC80');
plot(1:n+1, LB, '--', 'LineWidth', 2,'color','#EDB190');
plot(1:n+1, UB, '--','LineWidth', 2,'color','#EDB190');
plot(1:n+1, X,'LineWidth', 2,'color','b');
hold off
xlabel('Generation')
ylabel('Expectation')
legend('Estimation','','','X')
title('Confidence Interval')