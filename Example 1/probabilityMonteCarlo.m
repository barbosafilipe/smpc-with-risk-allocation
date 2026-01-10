function [Xu,Xl,BonfS] = probabilityMonteCarlo(x,nsim,mu,sigma,gU,gL,upp,low,N,nx,nw,A,G)

Gmap = zeros((N)*nx, N*nw);  % preallocate the zeros
for i = 1:N
    for j = 1:i
        Gmap((nx*(i-1)+1):(nx*i), (nw*(j-1)+1):(nw*j)) = (A^(i-j))*G;
    end
end


nx = size(x,1);
Xu = []; % to store the realizations in each prediction step uppee bound
Xl = []; % to store the realizations in each prediction step lower bound
rng('shuffle');
w = mu+chol(sigma)*randn(size(x,2),nsim); % create a random matrix normally distributed
xw = repmat(x(:),1,nsim) + Gmap*w; % compute the realizations
aux = []; % to help find the probabilities with Bonferroni
for i=1:size(gU,2)
    j = (i-1)*nx+1;
    Xu = [Xu; gU(i) nnz(xw(j,:)>upp)/nsim];
    Xl = [Xl; gL(i) nnz(xw(j,:)<low)/nsim];
    aux = [aux;xw(j,:)];
end
BonfS = [(sum(gU)+sum(gL)) nnz((aux(1,:)>upp) | (aux(2,:)>upp) | (aux(3,:)>upp) | (aux(4,:)>upp) | (aux(5,:)>upp) | (aux(6,:)>upp)  | (aux(1,:)<low) | (aux(2,:)<low) | (aux(3,:)<low) | (aux(4,:)<low) | (aux(5,:)<low) | (aux(6,:)<low))/nsim];
end