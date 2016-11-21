function [scores,resnorm] = ILSQ_Match(dA,U1,num_prediction,rr,cc,mask,integer)         
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
l0 = rand(cc,1);
opts = optimset('MaxFunEval',inf,'MaxIter',Inf,'display','off');

dA1 = reshape(dA,rr^2,1);

if nargin<7    %nargin<7, solve relaxed continuous lsq
    if nargin<6     %no mask
        scores = lsqlin(U1,dA1,[],[],[],[],zeros(cc,1),ones(cc,1),l0,opts);    %solve the relaxed continuous linear least squared problem
    else            %using mask
        mask = mask(:);
        [scores,resnorm] = lsqlin(bsxfun(@times,U1,mask),dA1.*mask,[],[],[],[],zeros(cc,1),ones(cc,1),l0,opts);    %solve the relaxed continuous linear least squared problem
        %[scores,resnorm] = cplexlsqlin(bsxfun(@times,U1,mask),dA1.*mask,[],[],[],[],zeros(cc,1),ones(cc,1),l0);    %solve using cplex, more stable (need cplex installation)
        %[scores,resnorm] = cplexlsqlin(bsxfun(@times,U1,mask),dA1.*mask,ones(1,cc),num_prediction,[],[],zeros(cc,1),ones(cc,1),l0);  %extra constraint sum_lambda < testnumber, not as good as no such constraint
    end
else            %nargin=7, solve integer lsq using CPLEX -- not feasible, too slow for variables > 1000
    %l = cplexlsqbilp(U1,dA1,ones(1,size(U1,2)),select);    %use IBM CPLEX to solve the integer least square problem, using selection number <= testnumber constraint
    scores = cplexlsqbilp(bsxfun(@times,U1,mask),dA1.*mask,ones(1,size(U1,2)),num_prediction);    %use mask
end
