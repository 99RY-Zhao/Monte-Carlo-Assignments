function [o] = tmarginal(lambda,theta,psi,tau,t,n,s,debug) 
    deltas = zeros(size(lambda));
    for i=1:length(lambda)
        deltas(i) = t(i+1)-t(i);
    end

    %============================================================%
    % DOES NOT WORK DUE TO OVER AND UNDERFLOWS                   %
    % We fixed this by using the exponential/logarithmic version %
    %============================================================%
    % o = gampdf(deltas,ones(size(lambda)),1./lambda);           %
    % o = o .* (lambda .^ (n-1));                                %
    % o = prod(o);                                               %
    %============================================================%
    o = exp(sum(log(lambda).*n + log(deltas) - lambda.*deltas));
    
    if debug
        fprintf("[Tm] %s: %f\n",s,o);
    end
end
