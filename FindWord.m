% Search for missing words (hyperlinks on alphabet).

s = 5000;
u = 2000;
nmissing = 2000;

rng('default');
%% load 10000 most frequent words
dataname = 'google-10000-english-usa.txt';
S = zeros(26, s); % incidence matrix of words
fid = fopen(['data/', dataname]);
tline = fgetl(fid);
i = 1;
while ischar(tline) && i<=s
    disp(tline)
    word = uint8(tline) - 96;
    for j = 1:length(word)
        S(word(j), i) = S(word(j), i) + 1;
    end
    tline = fgetl(fid);
    i = i + 1;
end
fclose(fid);


%% do statistics
lens = tabulate(sum(S, 1));
lens = lens(:, 3)/100;
abcs = sum(S, 2);
len_dist = makedist('Multinomial','probabilities',lens);
abc_dist = makedist('Multinomial','probabilities',abcs/sum(abcs));


%% generate 10000 random words
U = zeros(26, u);
for i = 1:u
    len = random(len_dist);
    word = random(abc_dist, 1, len);
    for j = 1:len
        U(word(j), i) = U(word(j), i) + 1;
    end
end


%% split training and testing data
perm = randperm(s);
dS = S(:, perm(1:nmissing));
S(:, perm(1:nmissing)) = [];
U = [U, dS];
A = S*S';
Ulabels = [zeros(1, u), ones(1, nmissing)];
%MATBoost(A, [], [], 1, U, nmissing, Ulabels);
%[Lambda, scores] = HLpredict(S, U, nmissing, 'SPHC', 1);
%nnz(Lambda'.*Ulabels)

% train a neural net
for runthis = 1:1
S2 = zeros(26, u);
for i = 1:u
    len = random(len_dist);
    word = random(abc_dist, 1, len);
    for j = 1:len
        S2(word(j), i) = S2(word(j), i) + 1;
    end
end
S = [S, S2];
ss = size(S, 2);
layers = [imageInputLayer([26 1 1], 'Normalization','none')
fullyConnectedLayer(32)
reluLayer
fullyConnectedLayer(32)
reluLayer
fullyConnectedLayer(2)
softmaxLayer
classificationLayer];
opts = trainingOptions('sgdm', 'InitialLearnRate', 0.1, 'MaxEpochs', 200, 'MiniBatchSize', 128, ...
    'LearnRateSchedule','piecewise', 'LearnRateDropFactor', 0.9, 'L2Regularization', 0, ...
    'ExecutionEnvironment', 'cpu');
net = trainNetwork(reshape(S, 26, 1, 1, ss), categorical([ones(1, s), zeros(1, u)]), layers, opts);
[~, scores] = classify(net, reshape(U, 26, 1, 1, u));
scores(:, 1) = [];
Lambda = zeros(size(scores, 1), 1);
[~,I] = sort(scores,1,'descend');
Lambda(I(1:nmissing)) = 1;    % only keep hl with top scores
nnz(Lambda'.*Ulabels)
end
