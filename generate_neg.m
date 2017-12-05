function neg = generate_neg(pos,mul)
%  Usage: for generating random negative examples with "mul" times size as pos.
%  --Input--
%  -pos: the positive example matrix you have, each column is an example
%  -mul: how many times more examples you want the generated neg to have
%  --Output--
%  -neg: the randomly generated negative examples
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
[row,col] = size(pos);

if nargin < 2
    mul = 1;
end
neg = zeros(row,col*mul);
for j = 1:mul
    for i = 1:col
        r = pos(:,i);
        num_meta = nnz(r);
        while 1
            c = zeros(row,1);
            perm = randperm(row);
            c(perm(1:num_meta)) = 1;
            %lia = ismember(c',pos','rows');    %check if c in US
            lia = 0;      %to save time
            if lia == 0
                i;
                break
            end
        end
        neg(:,i) = c;
    end
end

