function [newt,accepted] = mh(t,lambda,theta,psi,tau,n,f,rho,r,debug)
    accepted = length(t)-2;

    if r == "global"
        while true
            newt = rglob(t,rho);
            newn = getn(newt,tau,false);
            
            previous = newt(1);
            for i=1:(length(newt)-1)
                if newt(i+1) < previous
                    continue;
                end
            end

            % If we reach this code the proposal is valid.
            break;
        end

        % Evaluate transition probability
        pn = f(lambda,theta,psi,tau,newt,newn,"pn",debug(1));
        po = f(lambda,theta,psi,tau,t,n,"po",debug(2));
        prob = min(1,pn/po);
    
        if rand > prob
            % Negative! Don't perform transition
            newt = t;
            accepted = 0;
        end
    elseif r == "individual"
        newt = t;
        for i=2:length(lambda)
            update = rindv(newt,i,rho);
            while update < newt(i-1) || update > newt(i+1)
                update = rindv(newt,i,rho);
            end

            % Evaluate the update and the new incident distribution
            newt(i) = update;
            newn    = getn(newt,tau,false);

            % Evaluate transition probability
            pu = f(lambda,theta,psi,tau,newt,newn,"pu",debug(1));
            po = f(lambda,theta,psi,tau,t,n,"po",debug(2));
            prob = min(1,pu/po);

            if debug(3)
                fprintf("[");fprintf("%g, ",t(1:end-1));
                fprintf("%g] => [",t(end));fprintf("%g, ",newt(1:end-1));
                fprintf("%g]\n",newt(end));
            end

            if rand > prob
                % Negative! Don't perform transition
                newt(i) = t(i);
                accepted = accepted-1;
            end
        end
    else
        error("Unknown mode " + r);
    end
end

