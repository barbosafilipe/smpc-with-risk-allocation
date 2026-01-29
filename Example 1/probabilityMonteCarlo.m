function [Xu, Xl, BonfS] = probabilityMonteCarlo(x, nsim, mu, sigma, gU, gL, upp, low, N, nx, nw, A, G)

    % Preallocate Gmap with zeros
    Gmap = zeros(N*nx, N*nw);  
    for i = 1:N
        for j = 1:i
            Gmap((nx*(i-1)+1):(nx*i), (nw*(j-1)+1):(nw*j)) = (A^(i-j))*G;
        end
    end

    % Initialize variables
    nx = size(x, 1);
    Xu = []; % Store realizations for upper bound
    Xl = []; % Store realizations for lower bound
    rng('shuffle'); % Seed the random number generator
    w = mu + chol(sigma)*randn(size(x, 2), nsim); % Generate random matrix
    xw = repmat(x(:), 1, nsim) + Gmap*w; % Compute realizations
    aux = []; % Auxiliary variable for probabilities with Bonferroni

    % Calculate probabilities for upper and lower bounds
    for i = 1:size(gU, 2)
        j = (i-1)*nx + 1;
        Xu = [Xu; gU(i) nnz(xw(j,:) > upp) / nsim];
        Xl = [Xl; gL(i) nnz(xw(j,:) < low) / nsim];
        aux = [aux; xw(j,:)];
    end

    % Calculate Bonferroni statistic
    BonfS = [(sum(gU) + sum(gL)), ...
             nnz((aux(1,:) > upp) | (aux(2,:) > upp) | (aux(3,:) > upp) | ...
                 (aux(4,:) > upp) | (aux(5,:) > upp) | (aux(6,:) > upp) | ...
                 (aux(1,:) < low) | (aux(2,:) < low) | (aux(3,:) < low) | ...
                 (aux(4,:) < low) | (aux(5,:) < low) | (aux(6,:) < low))/nsim];
end
