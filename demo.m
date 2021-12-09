clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
ds = {'Caltech101-7'};
dsPath = '.\datasets\';
metric = {'ACC','nmi','Purity','Fscore','Precision','Recall','AR','Entropy'};
anchor_rate=[1 2 3 4 5 6 7];
d_rate = [1 2 3 4 5 6 7];

for dsi =1:1:1
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
    k = length( unique(Y));
    n = length(Y);
    
    %%
    for ichor = 1:length(anchor_rate)
        for id = 1:length(d_rate)
            tic;
            [A,W,Z,iter,obj,alpha,label] = algo_qp(X,Y,d_rate(id)*k,anchor_rate(ichor)*k);
            res = Clustering8Measure(Y, label);
            timer(ichor,id)  = toc;
            resall{ichor,id} = res;
            objall{ichor,id} = obj;
        end
    end
end


