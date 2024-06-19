function [t] = gett(limits,d)
    t = zeros(1,d+1);
    
    % Set the limits
    t(1) = limits(1);
    t(end) = limits(2);

    % Compute interval width
    w = (limits(2)-limits(1))/d;
    
    % Compute the breakpoints.
    for i=2:d
        t(i) = t(i-1) + w;
    end
end

