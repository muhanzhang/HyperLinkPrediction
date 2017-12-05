function [Lambda,scores] = CMM(hltrain,hltest,num_prediction,k,test_labels)
%  Usage: Main program of the Coordinated Matrix Minimization (CMM) algorithm for hyperlink prediction.
%  --Input--
%  hltrain: observed hyperlink columns (training hyperlinks)
%  hltest: candidate hyperlink columns (testing hyperlinks)
%  k: the number of latent factors used in nonnegative matrix factorization
%  others are the same as in HLpredict.m
%  --Output--
%  Lambda: logical column indicator vector (Lambda(i)==1 indicates the ith column in hltest is selected)
%  scores: column vector containing prediction scores of the columns in hltest
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
addpath(genpath('software/symnmf-master'));

if k == 'cv'
    if size(hltrain, 1) > 1000  % if hypernetwork is large, set k to default 30 to save time
        k = 30; 
    else   % otherwise, use cross validation to select k
        % Note that cross validation is tricky for transductive learning, not necessarily useful
        K = [10,20,30];
        folds = 2;  % should not use too many folds, which will make each fold too small
        ind = crossvalind('Kfold', size(hltrain,2), folds);
        AUC = zeros(length(K),folds);
        NMATCH = zeros(length(K),folds);
        display 'Cross validation begins...'
        for i = 1:folds
            hltest1 = [hltest, hltrain(:, ind==i)];
            hltrain1 = hltrain;
            hltrain1(:, ind==i) = [];
            test_labels1 = [zeros(size(test_labels)), ones(1, nnz(ind==i))];  % generate labels for this fold
            for j = 1:length(K)
                k = K(j);
                evalc('[nmatch, auc] = optimize(hltrain1,hltest1,nnz(ind==i),k,test_labels1);');
                AUC(j,i) = auc 
                NMATCH(j,i) = nmatch
            end
        end
        aAUC = mean(AUC,2);
        aNMATCH = mean(NMATCH,2);
        [~,I] = max(aNMATCH);  % use match number as cv criterion
        k = K(I);
    end
end

[nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels);
Lambda = logical(Lambda);

end




function [nmatch, auc, scores, Lambda] = optimize(hltrain,hltest,num_prediction,k,test_labels)
%  Usage: Alternate EM optimization between Lambda and W
%  --Output--
%  nmatch: number of true positive hyperlinks in the predictions
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%

A = hltrain * hltrain';  % the observed adjacency matrix
% reshape the test hyperlinks
[rr,cc] = size(hltest);
U1 = [];  % matrix of vectorized uu^T
for i=1:cc
    u = hltest(:,i);
    u = u*u';
    u = sparse(u(:));
    U1 = [U1, u];
end

% Optimization Settings
% settings of linear least square
opts = optimset('MaxFunEval',inf,'MaxIter',Inf,'display','on','algorithm','interior-point');
% settings of symmetric nonnegative matrix factorization
params = {}  
params.maxiter = 100;
params.debug = 0;
% global EM settings
max_iter = 100;  % maximum number of EM iterations
res = 100000;  % the current objective value

% optimization begins
scores = zeros(cc,1);   %initialize scores
Nmatch = [];
for iter = 1:max_iter

    old_res = res;

    Ai = A + (hltest*diag(scores)*hltest');  % A + U\Lambda U^T at ith iteration
    
    % symnnmf, the M step
    if iter == 1
        [W,~,res0] = symnmf_newton(Ai,k,params);
    else
        params.Hinit = W;  % if not the first iteration, optimize W starting from the old W
        [W,~,res0] = symnmf_newton(Ai,k,params);
    end
    WWT = W*W';
    
    % least square, the E step
    dA = WWT-A;
    [scores,res] = lsqlin(U1,dA(:),[],[],[],[],zeros(cc,1),ones(cc,1),[],opts);

    % results of the current iteration
    iter  % show iteration number
    res  % show current res
    Lambda = zeros(cc,1);
    [~,I] = sort(scores,1,'descend');
    Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
    nmatch = nnz(Lambda'.*test_labels)  % show number of true positive hls
    Nmatch = [Nmatch, nmatch];

    % early stopping
    if abs(res-old_res) < 1e-4
        break
    end
end
Nmatch
[~,~,~,auc] = perfcurve(test_labels,scores,true);
end
