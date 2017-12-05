function [Lambda,scores,hlpredicted] = HLpredict(hltrain,hltest,num_prediction,method,ith_experiment,test_labels)
%  Usage: A wrapper for different hyperlink prediction methods.
%         Call [Lambda,scores] = HLpredict(hltrain,hltest,num_prediction) to run in default mode.
%  --Input--
%  -hltrain: existing training hyperlinks, each row is a vertex, each column is a hyperlink. 
%           hltrain(i,j)==1 iff vertex i is in hyperlink j.
%  -hltest: the candidate hyperlinks, same format as hltrain
%  -num_prediction: how many number of hyperlinks to select as predictions
%  -method: 'CMM'--the Coordinated Matrix Minimization algorithm, default. 
%
%  -ith_experiment: a positive integer number indicating you are running the ith experiment 
%           (for parallelly running multiple experiments purpose, avoid file reading conflict)
%  -test_labels: the labels of hltest, used for showing intermediate results within an algorithm
%  --Output--
%  -Lambda: logical column indicator vector, Lambda(i)==1 iff hltest(:,i) is predicted as positive (0 otherwise)
%  -scores: column vector of prediction scores of hltest, higher means more likely to be positive
%  -hlpredicted: the predicted hyperlinks, same format as hltrain
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin < 4
    method = 'CMM';
end
if nargin < 5
    ith_experiment = 1;
end

A = hltrain * hltrain';   % project the hyperlinks into the vertex space to form the adjacency matrix (with weights)

method
switch method
    case 'CMM'
        [Lambda,scores] = CMM(hltrain,hltest,num_prediction,'cv',test_labels); % pass test_labels to CMM in order to show intermediate results
    case 'FM' 
        [Lambda,scores] = FM(hltrain,hltest,8,ith_experiment,num_prediction); % use the default k=8 in FM, not very sensitive to k
    case 'CN'
        [Lambda,scores] = HCN(hltest,A,num_prediction);
    case 'Katz'
        [Lambda,scores] = HKatz(hltrain,hltest,A,num_prediction);
    case 'HPLSF'
        [Lambda,scores] = HPLSF(hltrain,hltest,A,8,ith_experiment,num_prediction);
    case 'SHC'
        [Lambda,scores] = SHC(hltrain,hltest,num_prediction);
    case 'LR'
        [Lambda,scores] = LR(hltrain,hltest,num_prediction);
    case 'NN'
        [Lambda,scores] = NN(hltrain,hltest,num_prediction);
    case 'BS'
        [Lambda,scores] = BS(hltrain,hltest,num_prediction);
    case 'DPP'
        [Lambda,scores] = DPP(hltrain,hltest,num_prediction,30,test_labels);
end
hlpredicted = hltest(:,Lambda);
if nargin == 6  % if test_labels are provided, calculate and show the auc score
    [~,~,~,auc] = perfcurve(test_labels,scores,true);
    auc
end
