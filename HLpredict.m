function [I,scores,hlpredicted] = HLpredict(hltrain,hltest,num_prediction,method,ith_experiment,test_labels)
%  Usage: for hyperlink predictions in hypernetworks, 
%         call [I,scores] = HLpredict(hltrain,hltest,num_prediction) in default mode.
%  --Input--
%  -hltrain: existing training hyperlinks, each row is a vertex, each column is a hyperlink. 
%           hltrain(i,j)==1 iff the jth hyperlink involves vertex i.
%  -hltest: the test hyperlinks to predict which of them are positive, same format as hltrain
%  -num_prediction: how many number of hyperlinks do you want to keep as positive. If you input
%           value < 1, this will be used as a threshold. hls with scores > threshold will be
%           predicted as positive.
%  -method: 'MATBoost'--the Matrix Boosting algorithm, default. 
%
%  -ith_experiment: a positive integer number indicating you are running the ith experiment 
%           (for parallelly running multiple experiments purpose, avoid file reading conflict)
%  -validr: used specifically for 'MATBoost' alg, validation reactions
%  --Output--
%  -I: logical column indicator vector, I(i)==1 iff hltest(:,i) is predicted as positive (0 otherwise)
%  -scores: the column vector of scores that are assigned to each test hyperlink, higher scores 
%           indicates higher probabilities to be positive hyperlinks
%  -hlpredicted: the predicted hyperlinks, same format as hltrain
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin < 5
    ith_experiment = 1;
end

A = hltrain * hltrain';   % project the hyperlinks into the vertex space to form the adjacency matrix (with weights)
%A = A - diag(diag(A));   % remove self adjacency
k = 8;

switch method
    case 'MATBoost'
        [I,scores] = MATBoost(hltrain,hltest,num_prediction,'cv',test_labels);
    case 'CM'
        [I,scores] = MATBoost2(A,[],k,ith_experiment,hltest,num_prediction,1);
    case 'Greedy'
        [I,scores] = GreedyMatch(A,[],k,ith_experiment,hltest,num_prediction);
    case 'FM' 
        [I,scores] = FM_Match(hltrain, hltest,ones(size(hltest,2),1),k,ith_experiment,num_prediction);
    case 'HCN'
        [I,scores] = HCN(A,hltest,num_prediction);
    case 'HKatz'
        [I,scores] = HKatz(A,hltrain,hltest,num_prediction);
    case 'Submodular'
        [I,scores] = SubmodularMatch(A,[],k,ith_experiment,hltest,num_prediction);
    case 'HPLSF'
        [I,scores] = HPLSF(A,k,ith_experiment,hltrain,hltest,num_prediction);
    case 'SPHC'
        [I,scores] = SPHC(hltrain,hltest,num_prediction);
    case 'MDA'
        [I,scores] = MDA(hltrain,hltest,num_prediction);
    case 'LR'
        [I,scores] = LR(hltrain,hltest,num_prediction);
    case 'NN'
        [I,scores] = NN(hltrain,hltest,num_prediction);
    case 'BS'  % bayesian set
        [I,scores] = BS(hltrain,hltest,num_prediction);
end
hlpredicted = hltest(:,I);
if nargin == 6
    [~,~,~,auc] = perfcurve(valmatch,scores,true);
    auc
end
