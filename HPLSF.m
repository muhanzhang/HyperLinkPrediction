function [Lambda,scores] = HPLSF(hltrain,hltest,A,k,ith_experiment,num_prediction)
%  Hyperlink Prediction using Latent Social Features, from paper "Hyperlink Prediction in Hypernetworks Using Latent Social Features".
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
addpath(genpath('software/liblinear-2.1/matlab'));
[Feature,trainlabels,Feature1] = feature_extraction(hltrain,hltest,A,k,ith_experiment);  % entropy feature extraction

X = [Feature;Feature1];
Z = zscore(X);
Z(:,1:size(Z,2)-rank(Z)) = [];
Feature = sparse(Z(1:size(Feature,1),:));
Feature1 = sparse(Z(size(Feature,1)+1:end,:));

%% use Liblinear to train on the extracted features
[~, optim_c] = evalc('liblinear_train(trainlabels, Feature, ''-s 0 -C -q'');');
model = liblinear_train(trainlabels, Feature, sprintf('-s 0 -c %d -q', optim_c(1)));
[~, ~, scores] = liblinear_predict(ones(size(Feature1,1), 1), Feature1, model, '-b 1 -q');
l1 = find(model.Label == 1);   % position of label 1 in model.Label
scores = scores(:, l1);

Lambda = zeros(size(hltest,2),1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
Lambda = logical(Lambda);


function [Feature,trainlabels,Feature1] = feature_extraction(hltrain,hltest,A,k,ith_experiment,labels)
%  Usage: for generating features used in HPLSF or other feature-based hyperlink prediction methods.
%  supports two kinds of feature extraction methods: FM and MDS (Multidimensional scaling, default)
%  also supports directly outputing concatenated features and its labels (need to provide labels of hltest)
%  --Input--
%  -labels: a vector of 0/1 labels of hltest
%  -others are the same as in HLpredict.m
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%

method = 2;  % 1: use FM features (not supported anymore, please use the method 2 which is better), 2: use MDS features
%% erratum: the HLPSF results in the paper were gotten using method 1, should change to method 2 for better results

write_to_file = 0;  % whether to also write the generated


%% feature extraction using FM model (the w, v)
if method == 1  % this feature-generation method is no longer supported in this version of code.
    evalc('FM(A,test,k,ith_experiment,1);');   % run FM_SGD to output w, v features
    w=dlmread('data/FMmodelw.txt');
    v=dlmread('data/FMmodelv.txt');
    w0=dlmread('data/FMmodelw0.txt');
    % for testing if model is read correctly (by comparing the out_1.txt)
    for runthis = 1:1
        k=1;
        g=2;
        a = v(:,k);
        c = v(:,g);
        b = sum((a+c).^2 - a.^2-c.^2)/2 + w0 + w(k) + w(g);
        b1 = a'*c + w0 + w(k) + w(g);
        assert(b-b1<0.01);
    end
    V = [w';v];    %to combine v and w into a single matrix
end

%% feature extraction by min_Z ||D - ZZ'||, D is the shortest distance matrix.
if method == 2
    A = spones(A);
    D = eye(size(A));
    dis = 10*ones(size(A));
    for i = 1:5
        D = spones(D*A);
        tmp = D*i;
        dis = min(dis,tmp);
    end
    dis(dis==10) = 0;
    [V,D] = eigs(dis);
    k = min(k,size(D,1));
    V(:,k+1:end) = [];
    D = D(1:k,1:k);
    Z = V*(D^2);
    V = Z';
end
    
%% generate training examples
if write_to_file == 1
    fid = fopen(strcat('data/LatentFactorTrainingExamples.txt'),'w+');   % also directly output feature concatenations, potentially for variable-feature-length classifiers
end

n = 2;   %n times negative samples
trainlabels = [ones(1,size(hltrain,2)),zeros(1,n*size(hltrain,2))]';

% to generate positive training examples from trainr
Feature = [];
for i = 1:size(hltrain,2)
    r = hltrain(:,i);
    example = V(:,logical(r));
    feature = [];
    for j=1:size(V,1)
        feature = [feature,entropy(example(j,:))];
    end
    Feature = [Feature;feature];
    % output feature concatenations
    example = [1,reshape(example,1,[])];
    if write_to_file == 1
        fprintf(fid,strcat(num2str(example),'\r\n'));
    end
end
% to sample randomly from empties in hltest to form negative examples
for nn = 1:n
    for i = 1:size(hltrain,2)
        r = hltrain(:,i);
        num_meta = sum(r);
        
        while 1
            c = zeros(size(hltrain,1),1);
            perm = randperm(size(hltrain,1));
            c(perm(1:num_meta)) = 1;
            lia = 0;   %faster, not too much influence because sparsity
            %lia = ismember(c',hltrain','rows');    %check if c in US
            if lia == 0
                i;
                break
            end
        end
        example = V(:,logical(c));
        
        feature = [];
        for j=1:size(V,1)
            feature = [feature,entropy(example(j,:))];
        end
        Feature = [Feature;feature];
        example = [0,reshape(example,1,[])];
        if write_to_file == 1
            fprintf(fid,strcat(num2str(example),'\r\n'));
        end
    end
end
if write_to_file == 1
    fclose(fid);
end

%%  generate test examples
if write_to_file == 1
    fid2 = fopen(strcat('data/LatentFactorTestingExamples.txt'),'w+');
end
Feature1 = [];
for i = 1:size(hltest,2)
    r = hltest(:,i);
    example = V(:,logical(r));
    feature = [];
    for j=1:size(V,1)
        feature = [feature,entropy(example(j,:))];
    end
    Feature1 = [Feature1;feature];
    if nargin < 7
        labels = ones(size(hltest,2)) * 2;   % if hltest labels not provided, use "2" as labels instead
    end
    example = [labels(i),reshape(example,1,[])];
    if write_to_file == 1
        fprintf(fid2,strcat(num2str(example),'\r\n'));
    end
end
if write_to_file == 1
    fclose(fid2);
end
