function [Lambda,scores] = SPHC(hltrain,hltest,num_prediction)
%  Spectral Hypergraph Clustering, in "Learning with hypergraphs: Clustering, classification, and embedding"
%  scores = (I - xi * Theta)^-1 * y
%
%  *author: Muhan Zhang, Washington University in St. Louis
%% transpose all incidence matrices and calculate Theta
hltrain = hltrain';
hltest = hltest';
y = [ones(size(hltrain,1),1);zeros(size(hltest,1),1)];
H = [hltrain;hltest];
Dv = diag(sum(H,2));
De = diag(sum(H,1));
tmp = inv(Dv).^(1/2);
Theta = tmp*H*inv(De)*H'*tmp;

%% cross validation to choose 'a'
a = 1;  % typlically a = 1 is the best. To save time, one may skip the cross validation.
for runthis = 1:0
    K = 5;
    ind = crossvalind('Kfold', size(hltrain,1), K);
    A = [0.01,0.1,0.5,0.99,1];
    V = zeros(length(A),K);
    for i = 1:K
        yi = y;
        yi(ind==i) = 0;
        for j = 1:length(A)
            a = A(j);
            f = pinv(eye(size(H,1)) - a*Theta)*yi;
            scores = [f(ind==i);f(size(hltrain,1)+1:end)];
            Lambda = zeros(size(scores,1),1);
            [~,I] = sort(scores,1,'descend');
            if num_prediction > 1
                Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
                Lambda = logical(Lambda);
            else
                Lambda(scores>num_prediction) = 1;
                Lambda = logical(Lambda);
            end
            match = sum(Lambda(1:nnz(ind==i)));
            V(j,i) = match;   % use match number as cv criterion
            
        end
    end
    aV = mean(V,2);
    [~,I] = max(aV);
    a = A(I);
end

f = pinv(eye(size(H,1)) - a*Theta)*y;
scores = f(end-size(hltest,1)+1:end);
Lambda = zeros(size(hltest,1),1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end