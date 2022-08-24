%% Implementation of the BraneMF
clear all;close all;
addpath code

%% load networks


A(:,:,1) = load('Org4932BraneMF_w8_16June2022.mat').experimental;
A(:,:,2) = load('Org4932BraneMF_w8_16June2022.mat').neighborhood;
A(:,:,3) = load('Org4932BraneMF_w8_16June2022.mat').fusion;
A(:,:,4) = load('Org4932BraneMF_w8_16June2022.mat').cooccurence;
A(:,:,5) = load('Org4932BraneMF_w8_16June2022.mat').coexpression;
A(:,:,6) = load('Org4932BraneMF_w8_16June2022.mat').database;



%% setting parameters
dim = 500; %size of embedding
alpha2 = 1; %weignting factor


N = size(A,1); % number of vertices
M = size(A,3); % number of layers

B = (A(:,:,1)+A(:,:,2)+A(:,:,3)+A(:,:,4)+A(:,:,5)+A(:,:,6)/6);

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
    [p,fval,exitflag,output] = lbfgs(@comeig_lbfgs_A,p,N,M,A,D,Q,alpha,beta, options);
    
    P = reshape(p',N,N);
    fun = @comeig2_lbfgs_A;
    % solve Q while fixing P and D
    [q,fval,exitflag,output] = lbfgs(@comeig2_lbfgs_A,q,N,M,A,D,P,alpha,beta,options);
    Q = reshape(q',N,N);
    
    fun = @comeig_lbfgs_A;
    
    % evaluate the objective function
    
    cost(i) = comeig_lbfgs_A(p,N,M,A,D,Q,alpha,beta);
    
    plot(i,cost(i),'.r')
    hold on, drawnow
    
    % stopping criterion
    if i > 1 && abs(cost(i)-cost(i-1)) < 10^(-5)
        break
    end
    
end

for dim = [128,256,512,1024]
    for alpha2 = [0,0.25,0.50,0.75,1]

        v = P(:,1:dim);
        d = D(1:dim,1:dim);
        emb = v*(d^(alpha2));
        filename = strcat('Org4932_BraneMF_w8_d',num2str(dim),'_alpha_',num2str(alpha2),'.txt');
        writematrix(emb,filename,'Delimiter','tab')

    end
end



