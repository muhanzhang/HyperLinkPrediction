function [Lambda,scores] = MATBoost(hltrain,hltest,num_prediction,k,test_labels)
%  Usage: the main program of the MATBoost algorithm for hyperlink prediction. We suggest to
%         call MATBoost using HLpredict.m instead of directly calling this script.
%  --Input--
%  A: the weighted train adjacency matrix
%  test: no use, for consistency
%  k: the number of latent factors used in matrix factorization
%  max_iter: the maximum number of iterations allowed
%  others are the same as in HLpredict.m
%  --Output--
%  Lambda: logical column indicator vector (Lambda(i)==1 indicates the ith column in hltest is selected)
%  scores: column vector containing prediction scores of the columns in hltest
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
addpath(genpath('symnmf-master'));

if k == 'cv'
    K = [10,20,30];
    folds = 2;
    ind = crossvalind('Kfold', size(hltrain,2), folds);
    AUC = zeros(length(K),folds);
    NMATCH = zeros(length(K),folds);
    for i = 1:folds
        hltest1 = [hltest, hltrain(:, ind==i)];
        hltrain1 = hltrain;
        hltrain1(:, ind==i) = [];
        test_labels1 = [zeros(size(test_labels)), ones(1, nnz(ind==i))];
        for j = 1:length(K)
            k = K(j);
            [nmatch, auc] = optimize(hltrain1,hltest1,nnz(ind==i),k,test_labels1);
            AUC(j,i) = auc;   % use match number as cv criterion
            NMATCH(j,i) = nmatch;
        end
    end
    aAUC = mean(AUC,2);
    aNMATCH = mean(NMATCH,2);
    [~,I] = max(aNMATCH);
    k = K(I);
end
aAUC
aNMATCH
k
[nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels);

Lambda = logical(Lambda);

end





function [nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels)

A = hltrain * hltrain';
% reshape the test hyperlinks
[rr,cc] = size(hltest);
U1 = sparse(zeros(rr^2,cc));
for i=1:cc
    u = hltest(:,i);
    u = u*u';
    u = sparse(u(:));
    U1(:,i) = u;
end

% optimization settings
scores = zeros(cc,1);   %initialize scores
opt = statset('Maxiter',1000);
params = {}
params.maxiter = 100;
params.debug = 0;
Nmatch = [];
max_iter = 100;
res = 10000;
for iter = 1:max_iter
    
    old_res = res;
    if mod(iter,10)==0
        iter
    end
    Ak = A + (hltest*diag(scores)*hltest');
    
    % symnnmf
    if iter == 1
        [H, ~, res0] = symnmf_newton(Ak, k,params);
    else
        params.Hinit = H;
        [H, ~, res0] = symnmf_newton(Ak, k,params);
    end
    wdA = H*H';
    res0
    
    % integer least square matching
    [scores, res] = ILSQ_Match(wdA-A,U1,num_prediction,rr,cc);
    res
    
    Lambda = zeros(cc,1);
    [~,I] = sort(scores,1,'descend');
    Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
    nmatch = nnz(Lambda'.*test_labels)
    Nmatch = [Nmatch, nmatch];

    if abs(res-old_res) < 1e-4
        break
    end
    
end
Nmatch
[~,~,~,auc] = perfcurve(test_labels,scores,true);
end