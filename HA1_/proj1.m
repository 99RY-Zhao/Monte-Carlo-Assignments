addpath("Functions");
format compact;



%% EXERCISE 2 - Wind Turbine Power Estimation using MC
clearvars;
close all;
clc;

% Load the powercurve of the wind turbine.
load powercurve_V164;

% PARAMETERS
M = 12;        % Which month to inspect in graphs.
N = 10000;    % How many samples to take for MC estimation.
I = [3.5,25]; % Interval in which plot the functions for importance sampling.

% Define the wind distribution parameters for every month.
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
lambda = [10.6   9.7   9.2   8.0   7.8   8.1   7.8   8.1   9.1   9.9   10.6  10.6];
k      = [ 2.0   2.0   2.0   1.9   1.9   1.9   1.9   1.9   2.0   1.9    2.0   2.0];

% Generate samples for wind data.
Xmc = zeros(N,length(k));
Xtr = zeros(N,length(k));
Xis = zeros(N,length(k));
Vap = zeros(N/2,length(k));
Van = zeros(N/2,length(k));

for i=1:length(k)
    Xmc(:,i) = wblrnd(lambda(i),k(i),N,1);

    % Truncated MC
    Fa = wblcdf(I(1),lambda(i),k(i));
    Fb = wblcdf(I(2),lambda(i),k(i));
    Xtr(:,i) = wblinv(Fa + (Fb-Fa)*rand(N,1),lambda(i),k(i));

    % Antithetic MC
    U  = rand(N/2,1);
    Fa = wblcdf(3.5,lambda(i),k(i));
    Fb = wblcdf(25,lambda(i),k(i));
    I1 = U;
    I2 = (1-U);

    Vap(:,i) = wblinv(I1,lambda(i),k(i));
    Van(:,i) = wblinv(I2,lambda(i),k(i));
end

% Estimate the target with the Standard MC method.
figure(1), hold on, legend show;
title("Monte Carlo Estimation");
[Tmc,Tmce] = mcplot(Xmc,P,"Standard MC","r",M);
fn = @(x) trobj(x,I,lambda,k,P);
[Ttr,Ttre] = mcplot(Xtr,fn,"Truncated MC","g",M);

fn = @(x,m) cvobj(x,m,lambda,k,P);
[Tcv,Tcve] = mcplot({Xmc,Xmc},fn,"Control Variate MC","m",M);

% Estimate the target with the Importance Sampling method.
% Start by defining the parameters the IS algorithm is going to find
m = zeros(length(k));
s = zeros(length(k));
a = zeros(length(k));
d = cell(1,length(k));

for i=1:length(k)
    [m(i),s(i),a(i),d{i}] = is(N,P,[0 30],lambda(i),k(i));
    Xis(:,i) = normrnd(m(i),s(i),N,1);
end

fn = @(x) isobj(x,m,s,lambda,k,P);
[Tis,Tise] = mcplot(Xis,fn,"Importance Sampling","b",M);

fn = @(x,m) (P(x{1})+P(x{2}))/2;
[Tav,Tave] = mcplot({Vap,Van},fn,"Antithetic Sampling","#02fce3",M);

% Calculate the functions to plot for the month being inspected.
wpdf = wblpdf(d{M},lambda(M),k(M));
npdf = a(M)*normpdf(d{M},m(M),s(M));
objf = P(d{M})';
maxx = max([max(npdf) max(objf .* wpdf)]);

figure(2);
t = tiledlayout(1,3,"TileSpacing","compact");
title(t,"Functions for " + months{M}), nexttile, hold on;

yyaxis left,  plot(d{M},wpdf);
yyaxis right, plot(d{M},objf);
legend({'Probability distribution','Objective function'})
nexttile, hold on;

yyaxis left,  plot(d{M},npdf,"b");
yyaxis right, plot(d{M},wpdf .* objf,"r");
supplot(d{M},npdf,0,maxx,"blue");
supplot(d{M},(wpdf .* objf),0,maxx,"red");

legend({'Instrumental distribution','Product function'});
nexttile, hold on;

plot(d{M},(wpdf .* objf) ./ npdf);
legend({'Final integrand'});

% Plot the expected power output (with error) for every month.
figure(3), hold on;
title("Expected power output for each month");
errorbar(reordercats(categorical(months),months),Tmc,Tmce);
errorbar(reordercats(categorical(months),months),Ttr,Ttre);
errorbar(reordercats(categorical(months),months),Tcv,Tcve);
errorbar(reordercats(categorical(months),months),Tis,Tise);
errorbar(reordercats(categorical(months),months),Tav,Tave);
legend(["Standard MC","Truncated MC","Control Variate", ...
        "Importance Sampling","Antithetic Sampling"]);

% Calculate the probability of power generation for every month
prob = wblcdf(25,lambda,k) - wblcdf(3.5,lambda,k);
figure(4), hold on, ylim([0 1]);
title("Availability factor for each month");
bar(reordercats(categorical(months),months),prob);

% Use standard MC to estimate the total power flowing through a turbine
r = 1.225;
d = 164;
m = 3;

% Total power in the wind passing the turbine.
Ptot = (r*pi*(d^2)/8)*gamma(1 + m./k).*(lambda.^m);

% We already estimated P(v) and we reuse that result
Cm = Tav  ./ Ptot;
Ce = Tave ./ Ptot;

figure(5), hold on;
title("Expected power coefficient for each month");
errorbar(reordercats(categorical(months),months),Cm,Ce);

figure(6), hold on;
title("Expected capacity factor for each month");
bar(reordercats(categorical(months),months),Tav / 9.5e6);



%% EXERCISE 3 - Two wind turbines
clearvars;
close all;
clc;

% Load the powercurve of the wind turbine.
load powercurve_V164;

N = 50;     % Number of points in the grid for both v1 and v2
I = [3.5 25]; % Interval spanned by the grid for both v1 and v2 

% Compute all points of the grid
[M1,M2] = meshgrid(linspace(3.5,25,N));

% Weibull PDF and CDF with the selected parameters
w = @(x) wblpdf(x,9.13,1.96);
W = @(x) wblcdf(x,9.13,1.96);

% Joint probability distribution
a = 0.638;
p = 3;
q = 1.5;

f = @(x,y) w(x) .* w(y) .* (1 + a .* ...
    (1 - W(x) .^ p) .^ (q - 1) .* (1 - W(y) .^ p) .^ (q - 1) .* ...
    ((W(x) .^ p) .* (1 + p*q) - 1) .* ((W(y) .^ p) .* (1 + p*q) - 1));

% Obtain and visualize the joint PDF.
F = f(M1,M2);

figure(1), hold on, xlabel("V1"), ylabel("V2"), zlabel("PDF");
mesh(M1,M2,F,"EdgeColor","green");

% Define objective functions
phi1 = @(v1,v2) reshape(P(v1(:)),size(v1));                    % P(v1)
phi2 = @(v1,v2) reshape(P(v2(:)),size(v2));                    % P(v2)
phi3 = @(v1,v2) reshape(P(v1(:)).*P(v2(:)),size(v2));          % P(v1)P(v2)
phi4 = @(v1,v2) reshape((P(v1(:))+P(v2(:))).^2 ,size(v2));     % (P(v1)+P(v2))^2
phi5 = @(v1,v2) reshape((P(v1(:))+P(v2(:))) > 9.5e6,size(v2)); % Prob(P(v1)+P(v2) > 9.5)
phi6 = @(v1,v2) reshape((P(v1(:))+P(v2(:))) < 9.5e6,size(v2)); % Prob(P(v1)+P(v2) < 9.5)

% Define instrumental distribution parameters
mean1 = [12    10]; sigma1 = [28 6; 6 28];
mean2 = [10    12]; sigma2 = [28 6; 6 28];
mean3 = [12    12]; sigma3 = [16 1; 1 16];
mean4 = [12    12]; sigma4 = [20 2; 2 20];
mean5 = [10    10]; sigma5 = [28 2; 2 26];
mean6 = [5.6  5.6]; sigma6 = [20 0; 0 20];

% Define instrumental distributions
g1 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean1,sigma1),size(v1));
g2 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean2,sigma2),size(v1));
g3 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean3,sigma3),size(v1));
g4 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean4,sigma4),size(v1));
g5 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean5,sigma5),size(v1));
g6 = @(v1,v2) reshape(mvnpdf([v1(:) v2(:)],mean6,sigma6),size(v1));

% Plot the the instrumental distributions together with the products
% of the objective functions and probability distribution functions.
figure(2), hold on;
sgtitle("Comparison of Standard MC integrand with instrumental PDF");

% Calculate multipliers to rescale the plots and make them more readable.
% In particular, scale the instrumental PDF to have the same height of the
% product of the objective function and probability distribution function.
m1 = max(max(phi1(M1,M2).*F))/max(max(g1(M1,M2)));
m2 = max(max(phi2(M1,M2).*F))/max(max(g2(M1,M2)));
m3 = max(max(phi3(M1,M2).*F))/max(max(g3(M1,M2)));
m4 = max(max(phi4(M1,M2).*F))/max(max(g4(M1,M2)));
m5 = max(max(phi5(M1,M2).*F))/max(max(g5(M1,M2)));
m6 = max(max(phi6(M1,M2).*F))/max(max(g6(M1,M2)));

subplot(2,3,1), hold on, mesh(M1,M2,phi1(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m1*g1(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P1");

subplot(2,3,2), hold on, mesh(M1,M2,phi2(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m2*g2(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P2");

subplot(2,3,3), hold on, mesh(M1,M2,phi3(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m3*g3(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P3");

subplot(2,3,4), hold on, mesh(M1,M2,phi4(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m4*g4(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P4");

subplot(2,3,5), hold on, mesh(M1,M2,phi5(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m5*g5(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P5");

subplot(2,3,6), hold on, mesh(M1,M2,phi6(M1,M2).*F,"EdgeColor","g");
mesh(M1,M2,m6*g6(M1,M2),"EdgeColor","r");
legend(["Product","Instrumental function"]),title("Plot for P6");

% Detect if any points defy the constraint phi(x)*f(x) > 0 => g(x) > 0
display("Prob #1: " + sum((g1(M1,M2) == 0).*((phi1(M1,M2).*F) > 0),"all"));
display("Prob #2: " + sum((g2(M1,M2) == 0).*((phi2(M1,M2).*F) > 0),"all"));
display("Prob #3: " + sum((g3(M1,M2) == 0).*((phi3(M1,M2).*F) > 0),"all"));
display("Prob #4: " + sum((g4(M1,M2) == 0).*((phi4(M1,M2).*F) > 0),"all"));
display("Prob #5: " + sum((g5(M1,M2) == 0).*((phi5(M1,M2).*F) > 0),"all"));
display("Prob #6: " + sum((g6(M1,M2) == 0).*((phi6(M1,M2).*F) > 0),"all"));

display(newline);

figure(3), hold on;
sgtitle("Final integrand (rescaled)");

% Again we rescale the plots to have them more or less comparable.
subplot(2,3,1), mesh(M1,M2,phi1(M1,M2).*F./(m1*g1(M1,M2)),"EdgeColor","b");
title("Plot for P1");

subplot(2,3,2), mesh(M1,M2,phi2(M1,M2).*F./(m2*g2(M1,M2)),"EdgeColor","b");
title("Plot for P2");

subplot(2,3,3), mesh(M1,M2,phi3(M1,M2).*F./(m3*g3(M1,M2)),"EdgeColor","b");
title("Plot for P3");

subplot(2,3,4), mesh(M1,M2,phi4(M1,M2).*F./(m4*g4(M1,M2)),"EdgeColor","b");
title("Plot for P4");

subplot(2,3,5), mesh(M1,M2,phi5(M1,M2).*F./(m5*g5(M1,M2)),"EdgeColor","b");
title("Plot for P5");

subplot(2,3,6), mesh(M1,M2,phi6(M1,M2).*F./(m6*g6(M1,M2)),"EdgeColor","b");
title("Plot for P6");

% Performs the stochastic draws.
X1 = mvnrnd(mean1,sigma1,N);
X2 = mvnrnd(mean2,sigma2,N);
X3 = mvnrnd(mean3,sigma3,N);
X4 = mvnrnd(mean4,sigma4,N);
X5 = mvnrnd(mean5,sigma5,N);
X6 = mvnrnd(mean6,sigma6,N);

% Calculate the Importance Sampled Monte Carlo Estimate
T1 = mean(phi1(X1(:,1),X1(:,2)).*f(X1(:,1),X1(:,2))./g1(X1(:,1),X1(:,2)));
T2 = mean(phi2(X2(:,1),X2(:,2)).*f(X2(:,1),X2(:,2))./g2(X2(:,1),X2(:,2)));
T3 = mean(phi3(X3(:,1),X3(:,2)).*f(X3(:,1),X3(:,2))./g3(X3(:,1),X3(:,2)));
T4 = mean(phi4(X4(:,1),X4(:,2)).*f(X4(:,1),X4(:,2))./g4(X4(:,1),X4(:,2)));
T5 = mean(phi5(X5(:,1),X5(:,2)).*f(X5(:,1),X5(:,2))./g5(X5(:,1),X5(:,2)));
T6 = mean(phi6(X6(:,1),X6(:,2)).*f(X6(:,1),X6(:,2))./g6(X6(:,1),X6(:,2)));

display("E{P(V1)} = " + T1);
display("E{P(V2)} = " + T2);
display("E{P(V1)+P(V2)} = E{P(V1)}+E{P(V2)} = " + (T1+T2));
display("E{P(V1)P(V2)} = " + T3);
display("E{(P(V1)+P(V2))^2} = " + T4);

display(newline);

display("C{P(V1),P(V2)} = E{P(V1)P(V2)} - E{P(V1)}E{P(V2)} = " + (T3-(T2*T1)));
display("V{P(V1)+P(V2)} = E{(P(V1)+P(V2))^2} - E{P(V1)+P(V2)}^2 = " + (T4-(T1+T2)^2));
display("D{P(V1)+P(V2)} = sqrt(V{P(V1)+P(V2)}) = " + sqrt(T4-(T1+T2)^2));

display(newline);

% Estimate the Confidence interval for the two probabilities
ci5 = 1.96*sqrt((T5-T5^2)/N);
ci6 = 1.96*sqrt((T6-T6^2)/N);

display("P{P(V1)+P(V2) > 9.5} = " + T5 + " +/- " + ci5);
display("P{P(V1)+P(V2) < 9.5} = " + T6 + " +/- " + ci6);
display("P{P(V1)+P(V2) = 9.5} = " + (1-T5-T6));



%% CONVENIENCE FUNCTIONS
% Wrapper of the objective function for truncated MC that uses the
% appropriate parameters for each month.
function [y] = trobj(X,I,lambda,k,P)
    y = zeros(size(X))';
    
    for i=1:length(X)
        Fa = wblcdf(I(1),lambda(i),k(i));
        Fb = wblcdf(I(2),lambda(i),k(i));

        y(i) = P(X(i))*(Fb-Fa);
    end
end

% Wrapper of the objective function for CV that uses
% the appropriate parameters for each month.
function [Z] = cvobj(data,m,lambda,k,P)
    X = data{1};
    Y = data{2};

    E = gamma(1 + (1/k(m)))*lambda(m);
    C = -cov(P(X),Y);
    B = C(2,1)/var(Y);

    Z = P(X) + B.*(Y-E);
end

% Wrapper of the objective function for IS that uses the appropriate
% parameters for the normal and Weibull distributions for each month.
function [y] = isobj(X,m,s,lambda,k,P)
    y = zeros(size(X))';
    
    for i=1:length(X)
        objf = P(X(i));
        wpdf = wblpdf(X(i),lambda(i),k(i));
        npdf = normpdf(X(i),m(i),s(i));

        y(i) = (objf*wpdf)/npdf;
    end
end