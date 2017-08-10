function [Lambda,scores] = MATBoost(A,test,k,ith_experiment,hltest,num_prediction,valmatch)
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
%mask = logical(diag(ones(size(A,1),1)));
[rr,cc] = size(hltest);
max_iter = 100;
k=30;

% reshape the test hyperlinks
U1 = [];
for i=1:cc
    u = hltest(:,i);
    u = u*u';  % do not remove self adjacency here, because self adjacency will be kept in wdA
    u = sparse(u(:));    
    U1 = [U1,u];
end

scores = zeros(cc,1);   %initialize scores
%scores = rand(cc,1);   %initialize scores

Lambda = zeros(cc,1);   %initialize scores
DS = zeros(1,max_iter);  %record the ||newscores - scores||
Scores = [];
maxmatch = 0;
minres = 1000000;
cnt = 0;
oldres = 100000;

opt = statset('Maxiter',1000);
params = {}
params.maxiter = 100;
params.debug = 0;
Nmatch = [];

for iter = 1:max_iter
    if mod(iter,10)==0
        iter
    end
    
    % bias correction step
    %Ak = A + (hltest*diag(scores)*hltest').*(~mask);
    %Ak = A + (hltest*diag(Lambda)*hltest').*(~mask);
    %Ak = A + (hltest*diag(Lambda)*hltest');
    Ak = A + (hltest*diag(scores)*hltest');
    
    % completion step
    %[~,~,wdA] = evalc('FM(Ak,test,k,ith_experiment);');
    
    % nnmf
%     if iter == 1
%         [W,H] = nnmf(Ak,k);
%     else
%         [W,H] = nnmf(Ak,k,'w0',W,'h0',H,'algorithm','mult','options',opt);
%     end
%     wdA = W*H;
%     res1 = norm(Ak-wdA,'fro')^2
    
    % symnnmf
    
    if iter == 1
        [H, ~, res1] = symnmf_newton(Ak, k,params);
    else
        params.Hinit = H;
        [H, ~, res1] = symnmf_newton(Ak, k,params);
    end
    wdA = H*H';
    res1
    
    
    % matching step (ILSQ or Lasso, default ILSQ which is faster)
    %newscores = ILSQ_Match(wdA,U1,num_prediction,rr,cc,mask);



    [scores, res] = ILSQ_Match(wdA-A,U1,num_prediction,rr,cc,scores);

    %newscores = Lasso_Match(wdA,U1,num_prediction,rr,cc,mask,0.1);
    
    % stop criterion
    Lambda = zeros(cc,1);
    [~,I] = sort(scores,1,'descend');
    Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
    %Lambda = logical(Lambda);
    nmatch = nnz(Lambda'.*valmatch)
    Nmatch = [Nmatch, nmatch];
%     
%     if nmatch > maxmatch        % if stops because of stop criterion, calculate the mean scores before this iteration
%         maxmatch = nmatch;
%         scores = newscores;
%     end
    res

    if iter > 1
        [B0,I0] = sort(scores0,'descend');
        [B,I] = sort(scores,'descend');
        C = setdiff(I0(1:num_prediction),I(1:num_prediction));
        size(C,1)
        %display 'pause'
    end
    H0 = H;
    wdA0 = wdA;
    scores0 = scores;
    nmatch0 = nmatch;
    
    %if abs(res-res1)
%     if abs(res-oldres) < 1e-4
%         cnt = cnt + 1;
%     else
%         cnt = 0;
%     end
%     oldres = res;
%   
%     if cnt == 50
%         break
%     end
    
    
%     if res<minres
%         minres = res;
%         scores = newscores;
%     end
end
%scores = mean(Scores(:,1:max(end-1,1)),2);  % output the average scores
%scores = Scores(:,end-1);  % output the last learned scores
Nmatch

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
