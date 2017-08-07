function [scores,resnorm] = ILSQ_Match(dA,U1,num_prediction,rr,cc,l0)
%  Usage: a program used in MATBoost for solving the (integer) least square optimization
%         (the matching step for MATBoost).
%  --Input--
%  -dA: the completed adjacency matrix (the output of the completion step)
%  -U1: the reshaped test hyperlinks
%  -rr, cc: size of the original test hyperlinks incidence matrix
%  -mask: the mask passed by MATBoost
%  -integer: to solve the integer least square if integer~=None, otherwise solve the relaxed
%            continuous linear least square problem (default).
%  --Output--
%  -scores: prediction scores
%  -resnorm: the final value of the objective optimization problem
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
if nargin < 6
    l0 = rand(cc,1);
end
%opts = optimset('MaxFunEval',inf,'MaxIter',Inf,'display','on','algorithm','trust-region-reflective');
opts = optimset('MaxFunEval',inf,'MaxIter',Inf,'display','on','algorithm','interior-point');

%dA1 = reshape(dA,rr^2,1);
dA1 = dA(:);


%scores = cplexlsqbilp(U1,dA1,[],[]);    %use mask
[scores,resnorm] = lsqlin(U1,dA1,[],[],[],[],zeros(cc,1),ones(cc,1),l0,opts);    %solve the relaxed continuous linear least squared problem
%[scores,resnorm] = lsqlin(U1,dA1,[],[],[],[],zeros(cc,1),ones(cc,1),[],opts);    %solve the relaxed continuous linear least squared problem

%mask = mask(:);
%[scores,resnorm] = lsqlin(bsxfun(@times,U1,mask),dA1.*mask,[],[],[],[],zeros(cc,1),ones(cc,1),l0,opts);    %solve the relaxed continuous linear least squared problem
%[scores,resnorm] = cplexlsqlin(bsxfun(@times,U1,mask),dA1.*mask,[],[],[],[],zeros(cc,1),ones(cc,1),l0);    %solve using cplex, more stable (need cplex installation)
%[scores,resnorm] = cplexlsqlin(bsxfun(@times,U1,mask),dA1.*mask,ones(1,cc),num_prediction,[],[],zeros(cc,1),ones(cc,1),l0);  %extra constraint sum_lambda < testnumber, not as good as no such constraint


%l = cplexlsqbilp(U1,dA1,ones(1,size(U1,2)),select);    %use IBM CPLEX to solve the integer least square problem, using selection number <= testnumber constraint
%scores = cplexlsqbilp(bsxfun(@times,U1,mask),dA1.*mask,ones(1,size(U1,2)),num_prediction);    %use mask
%scores = cplexlsqbilp(U1,dA1);    %use mask

