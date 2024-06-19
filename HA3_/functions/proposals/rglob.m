% Update all the breakpoints with the same deviation
function [t] = rglob(t,rho)
    deviation = rho*(2*rand - 1);

    % Leave alone the first and last breakpoints.
    n = length(t) - 2;
    i = [2,n-1];

    t(i) = t(i) + deviation * ones(n);
end