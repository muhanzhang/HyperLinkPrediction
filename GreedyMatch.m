function [Lambda,scores] = GreedyMatch(A,test,k,ith_experiment,hltest,num_prediction)
%  Usage: a greedy version of MATBoost, which completes A first, and select test hyperlinks
%         which have top overlappings with A.
%  --Input--
%  A: the weighted train adjacency matrix
%  test: no use, for consistency
%  k: the number of latent factors used in matrix factorization
%  others are the same as in HLpredict.m
%  --Output--
%  Lambda: logical column indicator vector (Lambda(i)==1 indicates the ith column in hltest is selected)
%  scores: column vector containing prediction scores of the columns in hltest
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
% the completion step, get unweighted dA instead of wdA
[~,dA,~] = evalc('FM(A,test,k,ith_experiment);');

% the greedy matching step
[~,cc] = size(hltest);
scores = zeros(cc,1);
for i=1:cc
    r = hltest(:,i);
    rA = r*r';
    rA = rA - diag(diag(rA));  % remove self adjacency
    scores(i) = (nnz(rA.*dA))/nnz(rA);
end

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end