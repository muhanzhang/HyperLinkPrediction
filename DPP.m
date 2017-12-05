function [Lambda,scores] = DPP(hltrain,hltest,num_prediction,k,test_labels)
%  Usage: Use determinantal point process (DPP) to model the probability of a hyperlink (Beta)
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
addpath(genpath('software/symnmf-master'));

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
    aNMATCH = mean(NMATCH,2)
    [~,I] = max(aNMATCH);
    k = K(I)
end

[nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels);

Lambda = logical(Lambda);

end





function [nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels)

A = hltrain * hltrain';

% optimization settings
params = {}
params.maxiter = 100;
params.debug = 0;
[H, ~, res0] = symnmf_newton(A, k, params);
K = H * H'; % the kernel

% reshape the test hyperlinks
[rr,cc] = size(hltest);
scores = zeros(cc, 1);
for i=1:cc
    u = logical(hltest(:,i));
    Ku = K(u, u);
    scores(i) = nnz(u)*(det(Ku));
end
Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'ascend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
nmatch = nnz(Lambda'.*test_labels)
        
[~,~,~,auc] = perfcurve(test_labels,scores,true);
end
