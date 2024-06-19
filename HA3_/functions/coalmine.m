function [history,accrate] = coalmine(T,d,psi,rho,r,N,verbose)
    limits = [1658 1980]; % Define time interval.
    
    % Initialize the partition to d equidistant breakpoints.
    t = gett(limits,d);
    n = getn(t,T,verbose);
    
    % Initialize the random variables.
    theta  = gamrnd(2,1/psi);
    lambda = gamrnd(2,1/theta,1,d);

    history = cell(N,3);  
    accrate = 0;

    for i=1:N
        % Draw theta and lambda according to their marginals.
        theta  = gettheta(lambda,theta,psi,T,t,n);
        lambda = getlambda(lambda,theta,psi,T,t,n);
    
        % Get the breakpoints using MH
        [t,a] = mh(t,lambda,theta,psi,T,n,@tmarginal,rho,r,[false false false]);
        n = getn(t,T,false);
    
        % Update acceptance rate
        accrate = accrate + a;

        % Save history
        history{i,1} = theta;
        history{i,2} = lambda;
        history{i,3} = t;
    end

    accrate = accrate/(N*(length(t)-2));
end

