function [r,c]=convert_FM(train,outname)

datapath = strcat(pwd,'/data/');
%t = triu(train) - diag(diag(train));
[r,c,v] = find(train);  %assymetric matrix factorziation
a = [r,c,v];
dlmwrite(strcat(datapath,outname),a,'delimiter',';');
