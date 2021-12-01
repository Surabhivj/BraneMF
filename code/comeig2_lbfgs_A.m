function [f,g] = comeig2_lbfgs_A(q,N,M,L,D,P,alpha,beta)
% compute derivative with respect to P^(-1), which is used in the L-BFGS
% algorithm (treat P^(-1) as independent Q)
Q = reshape(q',N,N);
f = 0;

m = 1:M;
k = m(randperm(length(m)));

for i = k
tmp = 1/2*norm((L(:,:,i) - P*D*Q),'fro')^2;
f = f + tmp;
end

f = f + 1/2*alpha*norm(P,'fro')^2 + 1/2*alpha*norm(Q,'fro')^2 + 1/2*beta*norm(P*Q-eye(N),'fro')^2;

if nargout > 1
    G = zeros(size(Q));
    for i = k
        %tmp = (-1)*(L(:,:,i) - P*D*Q)*P'*D;
        tmp = (-1)*(D*P')*(L(:,:,i) - P*D*Q);
        G = G + tmp;
    end
    
    %G = G + (1/2.*Q) + beta.*P'*(P*Q-eye(N)) + ((-1)*(D*P')*(X - P*D*Q));
    G = G + (1/2.*Q)+ beta.*P'*(P*Q-eye(N));
    g = G(:)';
end