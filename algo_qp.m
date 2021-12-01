function [A,W,Z,iter,obj,alpha,label] = algo_qp(X,Y,d,numanchor)
% m      : the number of anchor. the size of Z is m*n.
% X      : n*di

%% initialize
maxIter = 50 ; % the number of iterations
IterMax = 50;

m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);

W = cell(numview,1);            % di * d
A = zeros(d,m);         % d  * m
Z = zeros(m,numsample); % m  * n
Z(:,1:m) = eye(m);

for i = 1:numview
   di = size(X{i},2); 
   W{i} = zeros(di,d);
   X{i} = mapstd(X{i}',0,1); % turn into d*n
end

alpha = ones(1,numview)/numview;
flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;
    
    %% optimize W_i
    AZ = A*Z; % since each view share the same A and Z.
    parfor iv=1:numview
        C = X{iv}*AZ';      
        [U,~,V] = svd(C,'econ');
        W{iv} = U*V';
    end

    %% optimize A
    sumAlpha = 0;
    part1 = 0;
    for ia = 1:numview
        al2 = alpha(ia)^2;
        sumAlpha = sumAlpha + al2;
        part1 = part1 + al2 * W{ia}' * X{ia} * Z';
    end
    [Unew,~,Vnew] = svd(part1,'econ');
    A = Unew*Vnew';
    
    %% optimize Z
    B = zeros(numsample,m);
    for iv=1:numview
        B = B+alpha(iv)^2*X{iv}'*W{iv}*A;
    end
    G = zeros(numsample,m);
    G(1:m,:) = eye(m);
    [label,~,P,~,~,term2] = coclustering_bipartite_fast1( B, G, numclass, sumAlpha, IterMax);
    Z = P';
    res = Clustering8Measure(Y, label);

    %% optimize alpha
    M = zeros(numview,1);
    for iv = 1:numview
        M(iv) = norm( X{iv} - W{iv} * A * Z,'fro');
    end
    Mfra = M.^-1;
    Q = 1/sum(Mfra);
    alpha = Q*Mfra;

    term1 = 0;
    for iv = 1:numview
        term1 = term1 + alpha(iv)^2 * norm(X{iv} - W{iv} * A * Z,'fro')^2;
    end
    
    obj(iter) = term1+term2;
    
    if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4 || iter>maxIter || obj(iter) < 1e-10)
        flag = 0;
    end
end
         
         
    
