function [o] = gettheta(lambda,theta,psi,tau,t,n)
    d = length(lambda);
    o = gamrnd(2+(2*d),1/(psi + sum(lambda)));
end

