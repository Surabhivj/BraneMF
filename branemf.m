%% Implementation of the BraneMF
clear all;close all;
addpath code

%% load networks

load('data/neighborhood.mat')
load('data/fusion.mat')
load('data/cooccurence.mat')
load('data/coexpression.mat')
load('data/experimental.mat')
load('data/database.mat')

A = full(net1);
A(:,:,2) = full(net2);
A(:,:,3) = full(net3);
A(:,:,4) = full(net4);
A(:,:,5) = full(net5);
A(:,:,6) = full(net6);

%% setting parameters
dim = 500; %size of embedding
alpha2 = 1; %weignting factor
N = size(A,1); % number of vertices
M = size(A,3); % number of layers

%% compute random walk matrix (make sure your system is compatible to execute python codes)

for i = 1:M
    L(:,:,i) = pyrunfile("rw_mat.py", "m", A = A(:,:,i), w = 3);  %A = adjacency matrix, w = window size
end

N = size(L,1); % number of vertices
M = size(L,3); % number of layers
%% initialization

B = (L(:,:,1)+L(:,:,2)+L(:,:,3)+L(:,:,4)+L(:,:,5)+L(:,:,6)/6);

niter = 100;
options = optimset('Display','iter');

figure()

[P,D,Q] = svds(B,6400,'largest');

p = P(:)';
q = Q(:)';
d =D(:)';

alpha = 0.1;
beta = 1;

%% begin joint decomposition

for i = 1:niter
    fun = @comeig_lbfgs_A;
    
    % solve P while fixing Q and D
    [p,fval,exitflag,output] = lbfgs(@comeig_lbfgs_A,p,N,M,L,D,Q,alpha,beta, options);
    
    P = reshape(p',N,N);
    fun = @comeig2_lbfgs_A;
    % solve Q while fixing P and D
    [q,fval,exitflag,output] = lbfgs(@comeig2_lbfgs_A,q,N,M,L,D,P,alpha,beta,options);
    Q = reshape(q',N,N);
    
    fun = @comeig_lbfgs_A;
    
    % evaluate the objective function
    
    cost(i) = comeig_lbfgs_A(p,N,M,L,D,Q,alpha,beta);
    
    plot(i,cost(i),'.r')
    hold on, drawnow
    
    % stopping criterion
    if i > 1 && abs(cost(i)-cost(i-1)) < 10^(-5)
        break
    end
    
end

v = P(:,1:dim);
d = D(1:dim,1:dim);
emb = v*(d^(alpha2));
filename = strcat('yeast_branemf_w_a_',num2str(alpha),'_b_',num2str(beta),'_d_',num2str(dim),'_alpha_',num2str(alpha2),'.txt');
writematrix(emb,filename,'Delimiter','tab')
