%  A program for generating the figures in paper.
%
%% *author: Muhan Zhang, Washington University in St. Louis

Model = {'iJO1366','iAF1260b','iAF692','iHN637','iIT341','iAB_RBC_283'};
Method = ['CMM','BS','SHC','HPLSF','FM','Katz','CN'];
MethodNum = [1:1:7];
Color = {'r','b',[248,147,29]/255,'c','m','g',[20,68,106]/255, 'k'};
LineStyle = {'-','-.','--','--','-.',':','-',':'};
MarkerType = {'x','d','+','v','o','*','^','.'};
w = 4;   %line width
w1 = 20;   %title font
w2 = 12;   %marker size
w3 = 20;   %legend size

cd ../result/
close all;

AUC = zeros(length(MethodNum), length(Model));
for mo = 1:6
    model = Model{mo};
    R = {};
    for me = 1:length(MethodNum)
        R{me} = load(sprintf('%s_%d.mat',model,MethodNum(me)));
        AUC(me, mo) = R{me}.average_AUC(end);
    end
    missingNumber = R{1}.missingNumber;
    
    for runthis = 1:0
    F1 = figure(mo);    %match number under different testnumber
    figure('Color',[1 1 1]);
    box on;
    grid on;
    hold on;
    F = {};
    for me = 1:length(MethodNum)
        F{me} = errorbar(missingNumber,R{me}.average_match_num,R{me}.std_match_num,[LineStyle{me} MarkerType{me}],'Color',Color{me},'LineWidth',w,'MarkerSize',w2);
    end
    F{me+1} = plot(missingNumber,R{1}.average_guess_match_num,[LineStyle{me+1} MarkerType{me+1}],'Color',Color{me+1},'LineWidth',w,'MarkerSize',w2);
    
    h = legend([F{1},F{2},F{3},F{4},F{5},F{6},F{7},F{8}],'CMM','BS','SHC','HPLSF','FM','Katz','CN','Random','Location','NorthWest');
    set(h,'Fontsize',w3);
    xlabel('Number of missing reactions','Fontsize',w1)
    ylabel('Number of recovered reactions','Fontsize',w1)
    title(model,'Fontsize',w1,'Interpreter','none');
    set(gca,'linewidth',1,'fontsize',w1,'xtick',missingNumber);
    axis([missingNumber(1) missingNumber(end) -inf inf]);
    mkdir Figures/;
    cd Figures/;
    evalc(sprintf('export_fig RM_%s -pdf -transparent -c[Inf,25,Inf,10]',Model{mo}));  % need to install export_fig
    cd ..;
    end
    
    % print auc and std
    tmp = '';
    Model(mo)
    for me = 1:length(MethodNum)
        tmp = strcat(tmp, ' & ', num2str(R{me}.average_AUC(end),'%.4f'),'$\pm$',num2str(R{me}.std_AUC(end),'%.4f'));
    end
    tmp;  % used for generating Latex Table lines
    
    for runthis = 1:0  % whether to draw AUC figures
        F2 = figure(mo+length(Model));    %AUC
        box on;
        grid on;
        hold on;
        F = {};
        for me = 1:length(Model)
            F{me} = plot(missingNumber,R{me}.average_AUC,[LineStyle{me} MarkerType{me}],'Color',Color{me},'LineWidth',w,'MarkerSize',w2);
        end
        F{me+1} = plot(missingNumber,0.5.*ones(1,length(missingNumber)),[LineStyle{me+1} MarkerType{me+1}],'Color',Color{me+1},'LineWidth',w,'MarkerSize',w2);
        
        h = legend([F{1},F{2},F{3},F{4},F{5},F{6},F{7},F{8}],'CMM','BS','SHC','HPLSF','FM','Katz','CN','Random','Location','Best');
        axis([missingNumber(1) missingNumber(end) -inf 0.9]);
        set(h,'Fontsize',w3);
        xlabel('Number of missing reactions','Fontsize',w1)
        ylabel('AUC','Fontsize',w1)
        title(model,'Fontsize',w1,'Interpreter','none');
        set(gca,'linewidth',1,'fontsize',w1,'xtick',missingNumber);
        cd Figures/;
        evalc(sprintf('export_fig AUC_%s -pdf -transparent -c[Inf,25,Inf,10]',Model{mo}));
        cd ..;
    end
    
    
    
end
cd ../utils/;
