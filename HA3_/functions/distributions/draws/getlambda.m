function [o] = getlambda(lambda,theta,psi,tau,t,n)
    deltas = zeros(size(lambda));

    for i=1:length(lambda)
        deltas(i) = t(i+1)-t(i);
    end

    o = gamrnd(2+n,1./(theta+deltas));
end

