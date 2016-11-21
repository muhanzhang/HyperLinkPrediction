function [scores,resnorm] = Lasso_Match(dA,U1,num_prediction,rr,cc,mask,lam)         
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
dA1 = reshape(dA,rr^2,1);


mask = mask(:);
U1 = bsxfun(@times,U1,mask);
dA1 = dA1.*mask;

lb = [zeros(cc,1)];
ub = [ones(cc,1)];

function [f,g] = obj(x)
f = (U1*x-dA1)'*(U1*x-dA1) + lam*sum(abs(x));
g = (2*U1'*U1*x - 2*U1'*dA1) + lam*sign(x);
end

x0 = [rand(cc,1)];

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',20000,'MaxIterations',1000,'SpecifyObjectiveGradient',true);
%options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',20000,'MaxIterations',10,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');

result = fmincon(@obj,x0,[],[],[],[],lb,ub,[],options);
scores = result(1:cc);

[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
end
