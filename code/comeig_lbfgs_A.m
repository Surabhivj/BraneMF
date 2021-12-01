function [f,g] = comeig_lbfgs_A(p,N,M,L,D,Q,alpha,beta)
% compute derivative with respect to P, which is used in the L-BFGS
% algorithm
P = reshape(p',N,N);
f = 0;

m = 1:M;
k = m(randperm(length(m)));

for i = k
tmp = 1/2*(norm((L(:,:,i) - P*D*Q),'fro'))^2;
f = f + tmp;
end

f = f + 1/2*alpha*norm(P,'fro')^2 + 1/2*alpha*norm(Q,'fro')^2 + 1/2*beta*norm(P*Q-eye(N),'fro')^2;

if nargout > 1
    G = zeros(size(P));
    for i = k
        tmp = (-1).*(L(:,:,i) - P*D*Q)*Q'*D;
        G = G + tmp;
      
    end
    G = G + (1/2.*P)+ beta.*(P*Q-eye(N))*Q';
    %G = G + (1/2.*P) + beta.*(P*Q-eye(N))*Q' + ((-1).*(X - P*D*Q)*Q'*D);
    g = G(:)';
end