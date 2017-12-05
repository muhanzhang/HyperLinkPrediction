function [Lambda,scores] = HCN(hltest,A,num_prediction)
%  Usage: Hypernetwork Common Neighbors Index for hyperlink prediction
%  CN = A*A
%  HCN(1,2,...,m) = (CN(1,2) + CN(1,3) + ... + CN(m-1,m)) / (m(m-1)/2),   1,2,...,m denote the vertices forming a hyperlink
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
[~,cc] = size(hltest);
train = spones(A);
sim = train*train;
sim = sim - diag(diag(sim));
sim = triu(sim);
scores = zeros(cc,1);
for i = 1:cc
    r = hltest(:,i);
    rA = r*r';
    rA = rA - diag(diag(rA));
    scores(i) = sum(sum(rA.*sim))/(nnz(r)*(nnz(r)-1)/2+1); % +1 to avoid dividing by zero
end

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
Lambda = logical(Lambda);
