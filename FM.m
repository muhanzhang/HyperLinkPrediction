function [dA,wdA] = FM(train,test,k,ith_experiment,SGD)
%  Usage: for solving the matrix completion problem, 
%         call [~,wdA] = FM(train,test,k,ith_experiment)
%         to get the weighted completed matrix wdA
%  --Input--
%  -train: the incomplete matrix
%  -test: no use in hl prediction (just for code consistency)
%  -k: number of latent factors in matrix factorization
%  -SGD: if SGD = 1, will use adaptive sgd for optimization (which can 
%        output parameters w, v). Otherwise, will use MCMC (default).
%  --Output--
%  -dA: the unweighted completed matrix, no need to use this
%  -wdA: the weighted completed matrix, which is what we want
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin<3
    k = 10;
end

FMtrain = sprintf('FMtrain_exp%d',ith_experiment);
FMtest = sprintf('FMtest_exp%d',ith_experiment);
FMconvert(train,FMtrain);      % call FMconvert.m to convert train, all to libfm format
all = ones(size(train));         % to predict all the possible links' scores
RC = FMconvert(all,FMtest);

optmethod = 1;
if nargin==5
    optmethod = 2;
end

switch optmethod
    case 1
        % to use MCMC optimization
        cd data;
        cmd = sprintf('./libFM -task r -train %s.libfm -test %s.libfm -dim "1,1,%d" -out out_%d.txt -iter 100',FMtrain,FMtest,k,ith_experiment);
        system(cmd);    %run libFM
        pred = dlmread(sprintf('out_%d.txt',ith_experiment));    %load the output file of libFM
        delete(sprintf('%s.libfm',FMtrain));
        delete(sprintf('%s.libfm',FMtest));
        delete(FMtrain);
        delete(FMtest);
        delete(sprintf('out_%d.txt',ith_experiment));
        cd ..;
    case 2
        % to use adaptive sgd with separate validation set
        [r,c,v] = find(train);
        rperm = randperm(length(r));
        val = rperm(1:floor(0.1*length(r)));    % split the train and validation data
        tra = rperm(floor(0.1*length(r))+1:end);
        validation = zeros(size(train));        % reassemble the matrices
        validation(sub2ind(size(validation),r(val),c(val))) = v(val);
        train = zeros(size(train));
        train(sub2ind(size(train),r(tra),c(tra))) = v(tra);
        FMconvert(train,FMtrain);
        FMvalidation = sprintf('FMvalidation_exp%d',ith_experiment);
        FMconvert(validation,FMvalidation);
        cd data;
        cmd = sprintf('./libFM -task r -train %s.libfm -test %s.libfm -validation %s.libfm -dim "1,1,%d" -out out_%d.txt -iter 100 -method sgda -learn_rate 0.001 -init_stdev 0.01',FMtrain,FMtest,FMvalidation,k,ith_experiment);
        system(cmd);                            % run libFM
        pred = dlmread(sprintf('out_%d.txt',ith_experiment));    % load the output file of libFM
        delete(sprintf('out_%d.txt',ith_experiment));
        delete(sprintf('%s.libfm',FMtrain));
        delete(sprintf('%s.libfm',FMtest));
        delete(sprintf('%s.libfm',FMvalidation));
        delete(FMtrain);
        delete(FMtest);
        delete(FMvalidation);
        cd ..;
end

%%
sim = zeros(size(train));
sim(sub2ind(size(sim),RC(:,1),RC(:,2))) = pred;    % assign FM predictions to sim
dsim = diag(diag(sim));     % diagonal elements of sim
sim = sparse((sim + sim'- dsim));  % NOTE: (sim + sim' - dsim) keeps diagonal elements

[~,dA,wdA] = CalcAUC(train,test,sim);
end
