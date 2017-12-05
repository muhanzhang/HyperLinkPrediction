function [Lambda,scores] = NN(hltrain,hltest,num_prediction)
%  Neural network for hl prediction.
%  Directly use columns of incidence matrix as positive training examples, and 
%  generate random columns as negative examples to train a neural network.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
S = full(hltrain);
lens = tabulate(sum(S, 1));
lens = lens(:, 3)/100;
abcs = sum(S, 2);
len_dist = makedist('Multinomial','probabilities',lens);
abc_dist = makedist('Multinomial','probabilities',abcs/sum(abcs));


U = full(hltest);
[m,s] = size(S);
[~,u] = size(U);

S2 = generate_neg(S,1);

S = [S, S2];
ss = size(S, 2);
layers = [imageInputLayer([m 1 1], 'Normalization','none')
fullyConnectedLayer(32)
reluLayer
fullyConnectedLayer(32)
reluLayer
fullyConnectedLayer(2)
softmaxLayer
classificationLayer];
opts = trainingOptions('sgdm', 'InitialLearnRate', 0.1, 'MaxEpochs', 200, 'MiniBatchSize', 128, ...
    'LearnRateSchedule','piecewise', 'LearnRateDropFactor', 0.9, 'L2Regularization', 0, ...
    'ExecutionEnvironment', 'cpu');
net = trainNetwork(reshape(S, m, 1, 1, ss), categorical([ones(1, s), zeros(1, size(S2,2))]), layers, opts);
[~, scores] = classify(net, reshape(U, m, 1, 1, u));
scores(:, 1) = [];
Lambda = zeros(size(scores, 1), 1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
Lambda = logical(Lambda);
