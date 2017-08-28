function [Lambda,scores] = BS(hltrain,hltest,num_prediction)
%  Bayesian set
%
%  *author: Muhan Zhang, Washington University in St. Louis
%% 
[row,N] = size(hltrain);
[~,col2] = size(hltest);

% prior
c = 2;
m = [hltrain, hltest];
m = mean(m,2);
alpha = c*m;
beta = c*(1-m);
alpha_t = alpha + sum(hltrain, 2);
beta_t = alpha + N - sum(hltrain, 2);

q = log(alpha_t) - log(alpha) - log(beta_t) + log(beta);
scores = hltest' * q;


Lambda = zeros(col2,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
