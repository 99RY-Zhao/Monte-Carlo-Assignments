function [p] = tauprior(lambda,t,n)
    d = length(t) - 1;
    e = 0;
    p = 1;

    for i=1:d
        p = p * lambda(i)^(n(i));
        e = e + lambda(i)*(t(i+1)-t(i));
    end

    e = exp(-e);
    p = e * p;
end