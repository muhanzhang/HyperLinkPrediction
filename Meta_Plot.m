%  This program is for plotting the results of Meta_RecoverS.m.
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%
Model = {'iJO1366','iAF1260b','e_coli_core','RECON1','iAT_PLT_636','iAF692','iHN637','iIT341'};
Method = ['MATBoost','SHC','HPLSF','FM','HKatz','HCN'];
MethodNum = [4,9,8,3,6,5];
Color = {'r','b','c','m',[248,147,29]/255,'g',[20,68,106]/255};
LineStyle = {'-','-.','--','--','-.',':',':'};
MarkerType = {'x','d','+','v','o','*','.'};
w = 4;   %line width
w1 = 20;   %title font
w2 = 12;   %marker size
w3 = 20;   %legend size

cd data/result/
close all;
for mo = 1:8
    
    
    model = Model{mo};
    
    R = {};
    for me = 1:length(MethodNum)
        R{me} = load(sprintf('%s_%d.mat',model,MethodNum(me)));
    end
    Testnumber = R{1}.Testnumber;
    
    F1 = figure(mo);    %match number under different testnumber
    box on;
    grid on;
    hold on;
    
    F = {};
    for me = 1:length(MethodNum)
        F{me} = plot(Testnumber,R{me}.average_match_num,[LineStyle{me} MarkerType{me}],'Color',Color{me},'LineWidth',w,'MarkerSize',w2);
    end
    F{me+1} = plot(Testnumber,R{1}.average_guess_match_num,[LineStyle{me+1} MarkerType{me+1}],'Color',Color{me+1},'LineWidth',w,'MarkerSize',w2);
    
    h = legend([F{1},F{2},F{3},F{4},F{5},F{6},F{7}],'MATBoost','SHC','HPLSF','FM','HKatz','HCN','Random','Location','NorthWest');
    set(h,'Fontsize',w3);
    xlabel('Missing number','Fontsize',w1)
    ylabel('Recovered number','Fontsize',w1)
    title(model,'Fontsize',w1,'Interpreter','none');
    set(gca,'linewidth',1,'fontsize',w1,'xtick',Testnumber);
    axis([Testnumber(1) Testnumber(end) -inf inf]);
    cd Figures/;
    evalc(sprintf('export_fig RM_%s -pdf -transparent -c[Inf,25,Inf,10]',Model{mo}));
    cd ..;
    
    
    
    F2 = figure(mo+length(Model));    %AUC
    box on;
    grid on;
    hold on;
    F = {};
    for me = 1:length(MethodNum)
        F{me} = plot(Testnumber,R{me}.average_AUC,[LineStyle{me} MarkerType{me}],'Color',Color{me},'LineWidth',w,'MarkerSize',w2);
    end
    F{me+1} = plot(Testnumber,0.5.*ones(1,length(Testnumber)),[LineStyle{me+1} MarkerType{me+1}],'Color',Color{me+1},'LineWidth',w,'MarkerSize',w2);
    
    h = legend([F{1},F{2},F{3},F{4},F{5},F{6},F{7}],'MATBoost','SHC','HPLSF','FM','HKatz','HCN','Random','Location','Best');
    axis([Testnumber(1) Testnumber(end) -inf 0.9]);
    set(h,'Fontsize',w3);
    xlabel('Missing number','Fontsize',w1)
    ylabel('AUC','Fontsize',w1)
    title(model,'Fontsize',w1,'Interpreter','none');
    set(gca,'linewidth',1,'fontsize',w1,'xtick',Testnumber);
    cd Figures/;
    evalc(sprintf('export_fig AUC_%s -pdf -transparent -c[Inf,25,Inf,10]',Model{mo}));
    cd ..;
end
cd ../..;