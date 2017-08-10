function [Lambda,scores] = HPLSF(wtrain,test,k,ith_experiment,hltrain,hltest,num_prediction)
%  Hyperlink Prediction using Latent Social Features, in "Hyperlink Prediction in Hypernetworks Using Latent Social Features".
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
addpath(genpath('/project/tantra/mh/LinkPrediction/software/liblinear-2.1/matlab'));
[Feature,trainlabels,Feature1] = feature_extraction(hltrain,hltest,wtrain,test,k,ith_experiment);   % for using shortest distance matrix decompostion features as in the original paper

X = [Feature;Feature1];
Z = zscore(X);
Z(:,1:size(Z,2)-rank(Z)) = [];
Feature = sparse(Z(1:size(Feature,1),:));
Feature1 = sparse(Z(size(Feature,1)+1:end,:));

%% use mnrfit to train lr model
% lr = mnrfit(Feature,trainlabels+1);
% scores = mnrval(lr,Feature1);
% scores(:,1) = [];

%% use Liblinear for lr
[~, optim_c] = evalc('liblinear_train(trainlabels, Feature, ''-s 0 -C -q'');');
model = liblinear_train(trainlabels, Feature, sprintf('-s 0 -c %d -q', optim_c(1)));
[~, ~, scores] = liblinear_predict(ones(size(Feature1,1), 1), Feature1, model, '-b 1 -q');
l1 = find(model.Label == 1);   % position of label 1 in model.Label
scores = scores(:, l1);

Lambda = zeros(size(hltest,2),1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end


