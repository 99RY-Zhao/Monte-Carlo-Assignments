function [mu,sigma,a,x] = is(N,P,range,lambda,k)
    % Define the sampled interval and the number of samples.
    x = linspace(range(1),range(2),N*(range(2)-range(1)));
    f = wblpdf(x,lambda,k);
    y = P(x)' .* f;
    a = max(y);

    % Template function for a normal distribution function.
    g = @(x,mu,sigma) normpdf(x,mu,sigma);

    champ = -1;
    sigma = -1;
    mu    = -1;

    for m=10:0.05:15
        for s=3:0.05:6
            res = var(y ./ g(x,m,s));
            if ~isnan(res) && (champ < 0 || res < champ)
                champ = res;
                mu = m;
                sigma = s;
            end
        end
    end

    display("Found N(" + mu + "," + sigma + ")");
end
