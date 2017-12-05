%  Main program of hyperlink prediction experiments on metabolic networks
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
clear all;

% mo: the model number, here a model means a reconstructed metabolic network
for mo = 3:3
        mo
    switch mo
        case 1
            model = 'iJO1366';
            missingNumber = 50:50:400;
        case 2
            model = 'iAF1260b';
            missingNumber = 50:50:400;
        case 3
            model = 'iAF692';
            missingNumber = 25:25:200;
        case 4
            model = 'iHN637';
            missingNumber = 25:25:200;
        case 5
            model = 'iIT341';
            missingNumber = 25:25:200;
        case 6
            model = 'iAB_RBC_283';
            missingNumber = 25:25:200;
    end
    load(sprintf('data/%s.mat',model),'Model');
    S = Model.S;    % stoichiometric matrix S, which is a signed incidence matrix
    S = full(spones(S));   % make S binary, now is a standard incidence matrix
    S_names = Model.rxns;   % reaction names (hyperplinks)
    M = Model.mets;   % metabolite names (nodes)
    U = full(spones(Model.US)); % universal reactions for this model, containing all columns of S
    U_names = Model.unrnames; % universal reactions' names

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Experiment begins. Change experimental settings here. %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numOfExperiment = 1;  % number of independent experiments

    % Method--1: CMM  2: BS  3: SHC  4: HPLSF  5: FM  6: Katz  7: CN  (8: LR  9: NN  10: DPP, not included in paper)
    Method = [1] % modify this list to compare more methods

    missingNumber = [200]  % Overwrites the default 25:25:200 or 50:50:400. To experiment under different missing numbers, comment this line

    for md = 1:length(Method)
        method = Method(md)
        len_missingNumber = length(missingNumber);
        Num_of_matched_reactions = zeros(len_missingNumber,numOfExperiment);
        Num_of_selected_reactions = zeros(len_missingNumber,numOfExperiment);
        Average_guess_match_num = zeros(len_missingNumber,numOfExperiment);
        AUC_of_Reaction_Prediction = zeros(len_missingNumber,numOfExperiment);
        for change_missingNumber = 1:len_missingNumber
            missnumber = missingNumber(change_missingNumber)  % the number of missing hyperlinks in the current iteration
            rseed = 2732*(missnumber^2-1);

            %poolobj = parpool(feature('numcores')); % uncomment this line and change the "for" in next line to "parfor" in order to run in parallel
            for ith_experiment = 1:numOfExperiment
                ith_experiment
                
                rand('seed',rseed*(ith_experiment^2-1));   % to reproduce multiple exp results
                
                perm = randperm(size(S,2));    % random split all hyperlinks into train and test
                train = sparse(S(:,perm(missnumber+1:end)));    % train reactions
                train_names = S_names(perm(missnumber+1:end));
                missing_names = S_names(perm(1:missnumber));

                tmp = strcmp(repmat(train_names,1,size(U_names,2)),repmat(U_names,size(train_names,1),1));  % compare train with U
                check_unrepeat = find(sum(tmp)==0);    
                test = U(:,check_unrepeat);  % keep those U which are not train reactions as test reactions
                test_names = U_names(1,check_unrepeat);
                
                check_labels = strcmp(repmat(missing_names,1,size(test_names,2)),repmat(test_names,size(missing_names,1),1));
                labels = sum(check_labels)>0;  % labels of test hyperlinks, 1: positive, 0: negative
                assert(nnz(labels>1)==0)

                % METHOD
                % Lambda: the 1/0 indicator vector indicating if each test hyperlink is predicted positive
                % scores: the soft indicator vector recording the confidence score of each test hyperlink being positive
                switch method
                    case 1
                        [Lambda,scores] = HLpredict(train,test,missnumber,'CMM',ith_experiment,labels);
                    case 2
                        [Lambda,scores] = HLpredict(train,test,missnumber,'BS',ith_experiment,labels);
                    case 3
                        [Lambda,scores] = HLpredict(train,test,missnumber,'SHC',ith_experiment,labels);
                    case 4
                        [Lambda,scores] = HLpredict(train,test,missnumber,'HPLSF',ith_experiment,labels);
                    case 5
                        [Lambda,scores] = HLpredict(train,test,missnumber,'FM',ith_experiment,labels);
                    case 6
                        [Lambda,scores] = HLpredict(train,test,missnumber,'Katz',ith_experiment,labels);
                    case 7
                        [Lambda,scores] = HLpredict(train,test,missnumber,'CN',ith_experiment,labels);
                    case 8
                        [Lambda,scores] = HLpredict(train,test,missnumber,'LR',ith_experiment,labels);
                    case 9
                        [Lambda,scores] = HLpredict(train,test,missnumber,'NN',ith_experiment,labels);
                    case 10
                        [Lambda,scores] = HLpredict(train,test,missnumber,'DPP',ith_experiment,labels);
                end

                predictions = test_names(Lambda');
                match = strcmp(repmat(missing_names,1,size(predictions,2)),repmat(predictions,size(missing_names,1),1)); 
                num_of_matched_reactions = nnz(match) % see how many predictions are real missing reactions
                Num_of_matched_reactions(change_missingNumber,ith_experiment) = num_of_matched_reactions;

                assert(nnz(Lambda)==missnumber)

                num_of_selected_reactions = nnz(Lambda);
                Num_of_selected_reactions(change_missingNumber,ith_experiment) = num_of_selected_reactions;
                average_guess_match_num = missnumber*missnumber/size(Lambda,1);
                Average_guess_match_num(change_missingNumber,ith_experiment) = average_guess_match_num;
                
                for output_predS = 1:0     % whether to output the true positive reactions and the recovered S + dS matrix
                    matched_reactions = predictions(logical(sum(match)))    % output the matched reactions
                    predS = full([train,test(:,logical(Lambda'))]);    % the recovered S + dS matrix (binary)
                end
                
                % calculate AUC
                [~,~,~,auc] = perfcurve(labels,scores,true);
                AUC_of_Reaction_Prediction(change_missingNumber,ith_experiment) = auc;
                
            end
            if exist('poolobj')
                delete(poolobj)
            end
            
        end

        % calculate mean results
        Recall = bsxfun(@times,Num_of_matched_reactions,1./missingNumber');
        Precision = Num_of_matched_reactions./Num_of_selected_reactions;
        average_recall = mean(Recall,2);
        std_recall = std(Recall,0,2);
        average_precision = mean(Precision,2);
        std_precision = std(Precision,0,2);
        average_guess_match_num = mean(Average_guess_match_num,2)
        std_guess_match_num = std(Average_guess_match_num,0,2);
        average_AUC = mean(AUC_of_Reaction_Prediction,2)
        std_AUC = std(AUC_of_Reaction_Prediction,0,2);
        average_match_num = mean(Num_of_matched_reactions,2)
        std_match_num = std(Num_of_matched_reactions,0,2);
        %sound(sin(2*pi*25*(1:4000)/100));  % make a sound in the end

        %% save results
        save(sprintf('result/%s_%d.mat',model,method),'average_match_num', 'std_match_num', 'average_guess_match_num', 'std_guess_match_num', 'average_AUC', 'std_AUC','average_recall','std_recall','average_precision','std_precision','missingNumber');
    end
end
