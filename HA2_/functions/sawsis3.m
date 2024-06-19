%% SIS IMPLEMENTATION - TAYLORED TO SELF AVOIDING RANDOM WALK PROBLEM
%  For convenience a local map is created and z,g are abstracted away
function [X,Y,Z,W] = sawsis3(n,k,g,resampling)
    if g == "naive"
        g = @g_naive;
    elseif g == "smart"
        g = @g_smart;
    else
        error("Wrong mode.");
    end

    X = zeros(k,n,"int16");
    Y = zeros(k,n,"int16");
    Z = zeros(k,n,"int16");
    W = zeros(k,n);

    % Initialize weights.
    W(1,:) = 1;

    % Store estimation of Cn for SISR.
    Cn = 1;

    for step=2:k
        for particle=1:n
            % If we run SISR this condition never evaluates true,
            if W(step-1,particle) == 0
                continue;
            end

            [Xnew,Ynew,Znew,Pnew] = g(X,Y,Z,step,particle,k);

            % No point in continuing if we crashed.
            fail = nnz((X(1:step-1,particle) == Xnew).*(Y(1:step-1,particle) == Ynew).*(Z(1:step-1,particle) == Znew));
            if fail == 0
                % Update the values
                X(step,particle) = Xnew;
                Y(step,particle) = Ynew;
                Z(step,particle) = Znew;
                
                % Update the weights
                W(step,particle) = W(step-1,particle)/Pnew;
            end
        end

        if resampling
            % Compute Cn
            Cn = Cn * mean(W(step,:));

            % Perform the resampling
            indexes = randsample(n,n,true,W(step,:));

            X(1:step,:) = X(1:step,indexes);
            Y(1:step,:) = Y(1:step,indexes);
            Z(1:step,:) = Z(1:step,indexes);

            % Reset the weights
            W(step,:) = ones(1,n);
        end
    end

    if resampling
        % We already computed Cn, we want to return that instead of W since
        % W is pretty useless in SISR since it's going to be ones(k,N)
        W = Cn;
    end
end



%% CONVENIENCE DISTRIBUTIONS
function [Xnew,Ynew,Znew,Pnew] = g_naive(X,Y,Z,step,particle,k)
    Xold = X(step-1,particle);
    Yold = Y(step-1,particle);
    Zold = Z(step-1,particle);

    choices = {[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]};
    action = randsample(choices,1);
    action = action{1};

    Xnew = Xold + action(1);
    Ynew = Yold + action(2);
    Znew = Zold + action(3);
    Pnew = 1/6;
end

function [Xnew,Ynew,Znew,Pnew] = g_smart(X,Y,Z,step,particle,k)
    Xnew = X(step-1,particle);
    Ynew = Y(step-1,particle);
    Znew = Z(step-1,particle);
    Pnew = 1;

    directions = {[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]};
    choices = {};

    for i=1:6
        choice = directions{i};

        % Calculate updated coordinates
        Xtest = Xnew + choice(1);
        Ytest = Ynew + choice(2);
        Ztest = Znew + choice(3);

        fail = nnz((X(1:step-1,particle) == Xtest).*(Y(1:step-1,particle) == Ytest).*(Z(1:step-1,particle) == Ztest));
        if fail == 0
            choices{end+1} = choice;
        end
    end

    if ~isempty(choices)
        action = randsample(choices,1);
        action = action{1};

        Xnew = Xnew + action(1);
        Ynew = Ynew + action(2);
        Znew = Znew + action(3);
        Pnew = Pnew/length(choices);
    end
end

