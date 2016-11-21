function [Lambda,scores] = LR(hltrain,hltest,num_prediction)
%  Logistic Regression for hl prediction, needs to install Liblinear package.
%  Directly use columns of incidence matrix as positive training examples, and 
%  generate random columns as negative examples to train a logitstic regression model.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%% 
[row,col] = size(hltrain);
[row,col2] = size(hltest);

%% generate negative training
neg = generate_neg(hltrain,1);
labels = [ones(col,1);zeros(size(neg,2),1)];
train = [hltrain, neg]';
train = sparse(train);

%% use liblinear for lr
[~, optim_c] = evalc('liblinear_train(labels, train, ''-s 0 -C -q'');');
model = liblinear_train(labels, train, sprintf('-s 0 -c %d -q', optim_c(1)));
[~, ~, scores] = liblinear_predict(ones(col2, 1), sparse(hltest'), model, '-b 1 -q');
l1 = find(model.Label == 1);   % position of label 1 in model.Label
scores = scores(:, l1);

Lambda = zeros(col2,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
