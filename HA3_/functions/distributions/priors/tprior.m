function [p] = tprior(t)
    d = length(t) - 1;
    p = 1;
    
    for i=1:d
        if t(i+1) < t(i)
            p = 0;
            return;
        end
        
        p = p * (t(i+1)-t(i));
    end
end