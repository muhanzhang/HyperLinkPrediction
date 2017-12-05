function [Lambda,scores] = FM(hltrain,hltest,k,ith_experiment,num_prediction)
%  Train a Factorization Machine on hyperlink columns.
%  Need executable "libFM" saved in data/FM_temp/
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
[row,col] = size(hltrain);
[row,col2] = size(hltest);
datapath = strcat(pwd,'/data/FM_temp/');       % path of the data folder
FMtrain = sprintf('FMtrain_exp%d',ith_experiment);

%% generate positive training examples
f1 = fopen(strcat(datapath,FMtrain),'w+');
for i = 1:col
    r = hltrain(:,i);
    ind = find(r);
    ind = [1,ind'-1];
    fprintf(f1,strcat(num2str(ind),'\r\n'));
end

%% generate negative training examples (give up, because performance is better if not generating negative)
Mul = 10;   % how many times more negative data
for runthis = 1:0
    for mul = 1:Mul
        for i = 1:col
            r = hltrain(:,i);
            num_meta = sum(r);
            while 1
                c = zeros(row,1);
                perm = randperm(size(hltrain,1));
                c(perm(1:num_meta)) = 1;
                %lia = ismember(c',hltest','rows');    % check if c in US
                lia = 0;      %to save time
                if lia == 0
                    i;
                    break
                end
            end
            ind = [0,find(c)'-1];
            fprintf(f1,strcat(num2str(ind),'\r\n'));
        end
    end
    fclose(f1);
end

%% generate testing
FMtest = sprintf('FMtest_exp%d',ith_experiment);
f2 = fopen(strcat(datapath,FMtest),'w+');
for i = 1:col2
    r = hltest(:,i);
    ind = find(r);
    ind = [1,ind'-1];
    fprintf(f2,strcat(num2str(ind),'\r\n'));
end
fclose(f2);

%% perform FM, need executable "libFM" saved in data/FM_temp/
cd data/FM_temp;
cmd = sprintf('python processFM.py %s',FMtrain);
system(cmd);
cmd = sprintf('python processFM.py %s',FMtest);
system(cmd);
cmd = sprintf('./libFM -task c -train %s.libfm -test %s.libfm -dim "1,1,%d" -out out_%d.txt -iter 1000',FMtrain,FMtest,k,ith_experiment);
evalc('system(cmd);');    % run libFM in silence mode
pred = load(sprintf('out_%d.txt',ith_experiment));    % load the predictions of libFM
delete(FMtrain);  % delete temporay files
delete(FMtest);
delete(sprintf('out_%d.txt',ith_experiment));
delete(sprintf('%s.libfm',FMtrain));
delete(sprintf('%s.libfm',FMtest));
cd ../..;

scores = pred;

Lambda = zeros(col2,1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:num_prediction)) = 1;    % only keep hl with top scores
Lambda = logical(Lambda);
