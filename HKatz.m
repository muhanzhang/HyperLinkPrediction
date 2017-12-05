function [Lambda,scores] = HKatz(hltrain,hltest,A,num_prediction)
%  Usage: Hypernetwork Katz Index for hyperlink prediction
%  Katz = (I - lambda * A)^-1 - I
%  HKatz(1,2,...,m) = (Katz(1,2) + Katz(1,3) + ... + Katz(m-1,m)) / (m(m-1)/2),   1,2,...,m denote the vertices forming a hyperlink
%
%  *author: Muhan Zhang, Washington University in St. Louis
%% cross validation to choose 'lambda'
lambda = 0.01;   %default value
for runthis = 1:1
    K = 5;
    ind = crossvalind('Kfold', size(hltrain,2), K);
    L = [0.001,0.005,0.01,0.1,0.5];
    V = zeros(length(L),K);
    for i = 1:K
        hl = hltrain(:,ind~=i);
        train = spones(hl*hl');
        train = train - diag(diag(train));
        hlt = [hltrain(:,ind==i),hltest];    % There are no negative hls in hltrain, maybe using hltest as neg is a better choice than random generation. 
        [~,cc] = size(hlt);
        for j = 1:length(L)
            lambda = L(j);
            sim = inv(sparse(eye(size(train,1))) - lambda * train);
            sim = sim - sparse(eye(size(train,1)));
            sim = triu(sim);
            scores = zeros(cc,1);
            for ii = 1:cc
                r = hlt(:,ii);
                rA = r*r';
                rA = rA - diag(diag(rA));
                scores(ii) = sum(sum(rA.*sim))/(nnz(r)*(nnz(r)-1)/2+1);
            end
            
            Lambda = zeros(cc,1);
            [~,I] = sort(scores,1,'descend');
            Lambda(I(1:num_prediction)) = 1;    % only keep hls with top scores
            Lambda = logical(Lambda);
            match = sum(Lambda(1:nnz(ind==i)));
            V(j,i) = match;    % cv criterion is match number
        end
    end
    V
    aV = mean(V,2);
    [~,I] = max(aV);
    lambda = L(I)
end

%% calculate HKatz using the best lambda
[~,cc] = size(hltest);
train = spones(A);
sim = inv(sparse(eye(size(train,1))) - lambda * train);
sim = sim - sparse(eye(size(train,1)));
sim = triu(sim);
scores = zeros(cc,1);
for i = 1:cc
    r = hltest(:,i);
    rA = r*r';
    rA = rA - diag(diag(rA));
    scores(i) = sum(sum(rA.*sim))/(nnz(r)*(nnz(r)-1)/2+1);
end

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
Lambda = logical(Lambda);
