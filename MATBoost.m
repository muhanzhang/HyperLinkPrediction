function [Lambda,scores] = MATBoost(A,test,k,ith_experiment,hltest,num_prediction,max_iter)
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
mask = spones(A)==0;
[rr,cc] = size(hltest);

% reshape the test hyperlinks
U1 = [];
for i=1:cc
    u = hltest(:,i);
    u = u*u';  % do not remove self adjacency here, because self adjacency will be kept in wdA
    u = sparse(u(:));    
    U1 = [U1,u];
end

scores = zeros(cc,1);   %initialize scores
DS = zeros(1,max_iter);  %record the ||newscores - scores||
Scores = [];

for iter = 1:max_iter
    
    % bias correction step
    Ak = A + (hltest*diag(scores)*hltest').*(~mask);
    
    % completion step
    [~,~,wdA] = evalc('FM(Ak,test,k,ith_experiment);');
    
    % matching step (ILSQ or Lasso, default ILSQ which is faster)
    newscores = ILSQ_Match(wdA,U1,num_prediction,rr,cc,mask);
    %newscores = Lasso_Match(wdA,U1,num_prediction,rr,cc,mask,0.1);
    
    % stop criterion
    ds = norm(newscores - scores);
    DS(iter) = ds;
    if ds>min(DS(1:max(iter-1,1)))        % if stops because of stop criterion, calculate the mean scores before this iteration
        break
    end
    scores = newscores;
    Scores = [Scores,scores];
end
scores = mean(Scores(:,1:max(end-1,1)),2);  % output the average scores
%scores = Scores(:,end-1);  % output the last learned scores

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
