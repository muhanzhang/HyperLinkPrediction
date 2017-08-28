% Search for missing dishes (hyperlinks on food species).
slCharacterEncoding('UTF-8')
u = 1000;
nmissing = 400;
dataname = 'chuancai.txt';
%dataname = 'yuecai.txt';

rng('default');
%% load data
% count number of different dishes
fid = fopen(['data/', dataname]);
tline = fgetl(fid);
i = 1;
j = 1;
dishmap = containers.Map;
foodmap = containers.Map;
while ischar(tline)
    disp(tline)
    words = split(tline, ':');
    dishname = words(1);
    dishname = dishname{1};
    if ~isKey(dishmap, dishname)
        dishmap(dishname) = i;
        i = i + 1;
    end
        
    foodnames = words(2);
    foodnames = foodnames{1};
    foodnames = split(foodnames, ',');    
    for k = 1:length(foodnames)
        if ~isKey(foodmap, foodnames(k))
            tmp = foodnames(k);
            tmp = tmp{1};
            foodmap(tmp) = j;
            j = j + 1;
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
N = foodmap.Count; % row number
M = dishmap.Count; % colunm number
S = zeros(N, M);
fid = fopen(['data/', dataname]);
tline = fgetl(fid);
while ischar(tline)
    words = split(tline, ':');
    dishname = words(1);
    dishname = dishname{1};
    foodnames = words(2);
    foodnames = foodnames{1};
    foodnames = split(foodnames, ',');    
    i = dishmap(dishname);
    for k = 1:length(foodnames)
        tmp = foodnames(k);
        tmp = tmp{1};
        j = foodmap(tmp);
        S(j,i) = 1;
    end
    tline = fgetl(fid);
end
fclose(fid);


%% do statistics
lens = tabulate(sum(S, 1));
lens = lens(:, 3)/100;
abcs = sum(S, 2);
len_dist = makedist('Multinomial','probabilities',lens);
abc_dist = makedist('Multinomial','probabilities',abcs/sum(abcs));


%% generate 10000 random words
U = zeros(N, u);
for i = 1:u
    while 1
    len = random(len_dist);
    word = random(abc_dist, 1, len);
    tmp = zeros(N,1);
    tmp(word) = 1;
    if ismember(tmp',S','rows') || ismember(tmp',U','rows')
        continue
    end
    for j = 1:len
        U(word(j), i) = U(word(j), i) + 1;
    end
    break
    end
end


%% split training and testing data
s = M;
for only_test_hls_more_than_one_node = 1:0
    [~,idx] = find(sum(S)==1);
    tmp = S(:,idx);
    S(:,idx) = [];
    s = size(S,2);
    S = [S,tmp];
end

perm = randperm(s);
dS = S(:, perm(1:nmissing));
S(:, perm(1:nmissing)) = [];
U = [U, dS];
Ulabels = [zeros(1, u), ones(1, nmissing)];

[~, scores] = MATBoost(S,U,nmissing,30,Ulabels);

%[Lambda, scores] = HLpredict(S, U, nmissing, 'BS', 1, Ulabels);
%nnz(Lambda'.*Ulabels)

nmissing^2/(u+nmissing)  % thoeretical random guess number


%% build reverse map and show created recipes
rev_foodmap = containers.Map('KeyType','double','ValueType','char');
allkeys = keys(foodmap);
for i = 1:length(allkeys)
    key = allkeys(i);
    key = key{1};
    rev_foodmap(foodmap(key)) = key;
end
[~, indices] = sort(scores, 1, 'descend');
for i = 1: length(indices)
    i
    idx = indices(i);
    rec = U(:, idx);
    if idx > u
        display 'This is a existing recipe.'
    else
        display 'This is a created recipe'
    end
    foods = find(rec==1);
    recipe = {};
    for j = 1:length(foods)
        recipe{j} = rev_foodmap(foods(j));
    end
    recipe
    waitforbuttonpress
end    

