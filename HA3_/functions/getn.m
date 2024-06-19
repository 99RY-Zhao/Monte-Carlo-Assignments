function [n] = getn(t,T,verbose)
    d = length(t)-1;
    n = zeros(1,d);

    for i=1:d
        a = T >= t(i);

        % Make sure to include the last year.
        if i < d
            b = T < t(i+1);
        else
            b = T <= t(i+1);
        end

        n(i) = sum(a .* b);
    end

    % Print some information about the subdivision.
    if verbose
        fprintf("We have the following breakpoints: %s\n" ,sprintf("%4.1f ",t));
        for i=1:d
            fprintf("  |-> Period %4.1f-%4.1f: %d\n",t(i),t(i+1),n(i));
        end
    end
end

