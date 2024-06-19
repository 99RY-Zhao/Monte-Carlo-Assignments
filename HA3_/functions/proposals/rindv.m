% Update every breakpoint with its own deviation
function [o] = rindv(t,i,rho)
    R = rho * (t(i+1)-t(i-1));
    deviation = R*(2*rand - 1);

    o = t(i) + deviation;
end