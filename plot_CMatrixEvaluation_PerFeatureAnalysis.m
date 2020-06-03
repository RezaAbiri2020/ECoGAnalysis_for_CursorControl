% This is the analysis for reviewers comments regarding to individual
% feature effect even before ltCLDA or on first day of ltCLDA:

%% Section 1: using the first day of ltCLDA to evaluate the C matrix for different
% features

clear all,
close all,
clc

%% Sestion 1-1: laoding the first day of ltCLDA
% experiment info:
% Long-Term CLDA in Center-Out and test performance in Center-Out 
expts = [];
expts(end+1).yymmdd = '20190521';
expts(end).CLDA_hhmmss = {'135008','141438'};
expts(end).perf_hhmmss = {'135943','141438'};

%% Section 1-2: Using CLDA blocks: calculating the new C matrices for individual feature and evaluate it through cross validation
AllFeatures=[];
AllKinematics=[];
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.CLDA_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            
              if length(TrialData.Events)==2 && TrialData.ErrorID==0% it means the cursor start from center and hit the target
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
              
                NeuralFeature=TrialData.NeuralFeatures(:,idx);
                AllFeatures=[AllFeatures, cell2mat(NeuralFeature)];
                AllKinematics=[AllKinematics, TrialData.CursorState(:,idx)];
                
              end       
        end   
    end
end
% performig regress for individual feature to calculate the C matrices

Day0_CLDA_Blocks_R2{1,1}='Features';
Day0_CLDA_Blocks_R2{2,1}='PhaseofDelta';
Day0_CLDA_Blocks_R2{3,1}='Delta';
Day0_CLDA_Blocks_R2{4,1}='Theta';
Day0_CLDA_Blocks_R2{5,1}='Alpha';
Day0_CLDA_Blocks_R2{6,1}='Beta';
Day0_CLDA_Blocks_R2{7,1}='LowGamma';
Day0_CLDA_Blocks_R2{8,1}='HighGamma';

Day0_CLDA_Blocks_R2{1,2}='Posx';
Day0_CLDA_Blocks_R2{1,3}='Posy';
Day0_CLDA_Blocks_R2{1,4}='Velx';
Day0_CLDA_Blocks_R2{1,5}='Vely';
Day0_CLDA_Blocks_R2{1,6}='Const';

for Fe=1:7
    Input_dynamics=AllFeatures(128*(Fe-1)+1:(128*(Fe-1)+128),:);
    
    for Kin=1:5
        Output_Dynamics=AllKinematics(Kin,:);
        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
        Day0_CLDA_Blocks_R2{1+Fe,1+Kin}=Rsq_Pre_output;
        
    end
end

%% Section 1-3: Using Performance blocks: calculating the new C matrices for individual feature and evaluate it through cross validation
AllFeatures=[];
AllKinematics=[];
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.perf_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.perf_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.perf_hhmmss{1,j},'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            
              if length(TrialData.Events)==2 && TrialData.ErrorID==0% it means the cursor start from center and hit the target
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
              
                NeuralFeature=TrialData.NeuralFeatures(:,idx);
                AllFeatures=[AllFeatures, cell2mat(NeuralFeature)];
                AllKinematics=[AllKinematics, TrialData.CursorState(:,idx)];
                
              end       
        end   
    end
end
% performig regress for individual feature to calculate the C matrices

Day0_Fixed_Blocks_R2{1,1}='Features';
Day0_Fixed_Blocks_R2{2,1}='PhaseofDelta';
Day0_Fixed_Blocks_R2{3,1}='Delta';
Day0_Fixed_Blocks_R2{4,1}='Theta';
Day0_Fixed_Blocks_R2{5,1}='Alpha';
Day0_Fixed_Blocks_R2{6,1}='Beta';
Day0_Fixed_Blocks_R2{7,1}='LowGamma';
Day0_Fixed_Blocks_R2{8,1}='HighGamma';

Day0_Fixed_Blocks_R2{1,2}='Posx';
Day0_Fixed_Blocks_R2{1,3}='Posy';
Day0_Fixed_Blocks_R2{1,4}='Velx';
Day0_Fixed_Blocks_R2{1,5}='Vely';
Day0_Fixed_Blocks_R2{1,6}='Const';

for Fe=1:7
    Input_dynamics=AllFeatures(128*(Fe-1)+1:(128*(Fe-1)+128),:);
    
    for Kin=1:5
        Output_Dynamics=AllKinematics(Kin,:);
        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
        Day0_Fixed_Blocks_R2{1+Fe,1+Kin}=Rsq_Pre_output;
        
    end
end

%% Section 2: using the imagined blocks to evaluate the C matrix for different
% features 
% Can estimation be better for C using multiple days of imagined blocks?

clear all,
close all,
clc

%% Sestion 2-1: laoding the imagined data
% experiment info:
% imagined in Center-Out 
expts = [];

expts(end+1).yymmdd = '20190501';
expts(end).Imagined_hhmmss = {'133745'};

expts(end+1).yymmdd = '20190507';
expts(end).Imagined_hhmmss = {'105331'};

expts(end+1).yymmdd = '20190510';
expts(end).Imagined_hhmmss = {'104601','132311'};

expts(end+1).yymmdd = '20190514';
expts(end).Imagined_hhmmss = {'111822'};

expts(end+1).yymmdd = '20190515';
expts(end).Imagined_hhmmss = {'105804','132722'};

%% Section 2-2: calculating the new C matrices for individual/all feature and evaluate it through cross validation
% also adding imagined blocks/days to calculation to see the if C got better

MultiSession_Features=[];
MultiSession_Kinematics=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
        
    for j=1:length(expt.Imagined_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.Imagined_hhmmss{1,j})
        
        % go through datafiles in imagined blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.Imagined_hhmmss{1,j},'Imagined');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        PerSession_Features=[];
        PerSession_Kinematics=[];
        
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            idx=find(TrialData.Time>=TrialData.Events(2).Time);
            NeuralFeature=TrialData.NeuralFeatures(:,idx);
            PerSession_Features=[PerSession_Features, cell2mat(NeuralFeature)];
            PerSession_Kinematics=[PerSession_Kinematics, TrialData.CursorState(:,idx)];
            
        end
        
        % performig regress for individual feature to calculate the C matrices
        for Fe=1:7
            Input_dynamics=PerSession_Features(128*(Fe-1)+1:(128*(Fe-1)+128),:);
            
            for Kin=1:5
                Output_Dynamics=PerSession_Kinematics(Kin,:);
                Figures=0;
                Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
                % the row is feature the column is the kin 
                R2Day(i).PerSession(j).Fe(Fe,Kin)=Rsq_Pre_output;
   
            end
        end
        
        %for performig regress for all features together per session to
        %calculate the C matrix
        Input_dynamics=PerSession_Features;
        for Kin=1:5
            Output_Dynamics=PerSession_Kinematics(Kin,:);
            Figures=0;
            Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
            % the row is feature the column is the kin
            R2Day(i).PerSession(j).AllFe(Kin)=Rsq_Pre_output;
                       
        end
        
        % for performig regress for individual feature to the current session to calculate the C matrices
        MultiSession_Features=[MultiSession_Features, PerSession_Features];
        MultiSession_Kinematics=[MultiSession_Kinematics, PerSession_Kinematics];
        
        for Fe=1:7
            Input_dynamics=MultiSession_Features(128*(Fe-1)+1:(128*(Fe-1)+128),:);
            
            for Kin=1:5
                Output_Dynamics=MultiSession_Kinematics(Kin,:);
                Figures=0;
                Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
                % the row is feature the column is the kin
                R2Day(i).MultiSession(j).Fe(Fe,Kin)=Rsq_Pre_output;
                
            end
        end
        
        %for performig regress for all features together to the current session
        Input_dynamics=MultiSession_Features;
        for Kin=1:5
            Output_Dynamics=MultiSession_Kinematics(Kin,:);
            Figures=0;
            Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
            % the row is feature the column is the kin
            R2Day(i).MultiSession(j).AllFe(Kin)=Rsq_Pre_output;
            
        end
        
    end
      
end

FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};
Kin_Parameter={'Posx','Posy','Velx','Vely'};

% plot the changes in R2 per feature across session; per session
for Fe=1:7
    
    figure(Fe);
    set(gcf, 'Position', [100, 100, 600, 800]);
    suptitle({[FeatureName{Fe},'; R2(in kinematics prediction) value per session']
        ['using data within each session']})
    
    Values=[];
    for Kin=1:4
        Counter=0;
        for i=1:5
            for j=1:length(R2Day(i).PerSession)
                Counter=Counter+1;
                Values(Kin,Counter)=R2Day(i).PerSession(j).Fe(Fe,Kin);           
            end  
        end
    end
    
    for Kin=1:4
        subplot(4,1,Kin)
        %title(Kin_Parameter{Kin})
        plot(Values(Kin,:),'ob-','linewidth',1)
        xlabel('Session')
        set(gca, 'fontsize',10)
        xlim([0 8])
        legend(Kin_Parameter{Kin})
    end
   
    %HighQualityFigs(['R2PerSession_',FeatureName{Fe}])
end


% plot the changes in R2 per feature across session; multi session
for Fe=1:7
    
    figure(Fe);
    set(gcf, 'Position', [100, 100, 600, 800]);
    suptitle({[FeatureName{Fe},'; R2(in kinematics prediction) value per session']
        ['using data up to current session']})
    
    Values=[];
    for Kin=1:4
        Counter=0;
        for i=1:5
            for j=1:length(R2Day(i).MultiSession)
                Counter=Counter+1;
                Values(Kin,Counter)=R2Day(i).MultiSession(j).Fe(Fe,Kin);           
            end  
        end
    end
    
    for Kin=1:4
        subplot(4,1,Kin)
        %title(Kin_Parameter{Kin})
        plot(Values(Kin,:),'ob-','linewidth',1)
        xlabel('Session')
        set(gca, 'fontsize',10)
        xlim([0 8])
        legend(Kin_Parameter{Kin})
    end
   
    %HighQualityFigs(['R2MultiSession_',FeatureName{Fe}])
end

% plot the changes in R2 for all features across session; per session
% calculation 
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);
suptitle({'All features; R2(in kinematics prediction) value per session for all features'
    'using data only per session'})

Values_PerSession=[];
for Kin=1:4
    Counter=0;
    for i=1:5
        for j=1:length(R2Day(i).PerSession)
            Counter=Counter+1;
            Values_PerSession(Kin,Counter)=R2Day(i).PerSession(j).AllFe(Kin);
        end
    end
end

for Kin=1:4
    subplot(4,1,Kin)
    %title(Kin_Parameter{Kin})
    plot(Values_PerSession(Kin,:),'ob-','linewidth',1)
    xlabel('Session')
    set(gca, 'fontsize',10)
    xlim([0 8])
    legend(Kin_Parameter{Kin})
end

%HighQualityFigs('R2PerSession_AllFeatures')

% plot the changes in R2 for all features across session; multi session
figure;
set(gcf, 'Position', [100, 100, 1000, 800]);
suptitle({'All features; R2(in kinematics prediction) value per session for all features'
    'using data up to current session'})

Values_MultiSession=[];
for Kin=1:4
    Counter=0;
    for i=1:5
        for j=1:length(R2Day(i).MultiSession)
            Counter=Counter+1;
            Values_MultiSession(Kin,Counter)=R2Day(i).MultiSession(j).AllFe(Kin);
        end
    end
end

for Kin=1:4
    subplot(4,1,Kin)
    %title(Kin_Parameter{Kin})
    plot(Values_MultiSession(Kin,:),'ob-','linewidth',1)
    xlabel('Session')
    set(gca, 'fontsize',10)
    xlim([0 8])
    legend(Kin_Parameter{Kin})
end

%HighQualityFigs('R2MultiSession_AllFeatures')

% for paper Figure 1
figure;
set(gcf, 'Position', [100, 100, 400, 110]);
plot(Values_PerSession(3,:),'ob--','linewidth',1,'MarkerEdgeColor','b')
hold on
plot(Values_PerSession(4,:),'or--','linewidth',1,'MarkerEdgeColor','r')
hold on
plot(Values_MultiSession(3,:),'ob-','linewidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
hold on
plot(Values_MultiSession(4,:),'or-','linewidth',1,'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
set(gca,'fontsize',10)
xticks([1:8])
xlim([.5, 7.5])
yticks([0:.2:.8])
ylim([0, .8])
box off



%% Section 3: using the last day of imagined blocks to evaluate the C matrix for different
% features


clear all,
close all,
clc

%% Sestion 3-1: laoding the last day of imagined blocks
% experiment info:
% CLDA in Center-Out and test performance in Center-Out 
expts = [];
expts(end+1).yymmdd = '20190515';
expts(end).Imagined_hhmmss = {'105804','132722'};
expts(end).CLDA_hhmmss = {'110735'};
expts(end).perf_hhmmss = {'112302','113447'};

%% Section 3-2: Using imagined blocks: calculating the new C matrices for individual feature and evaluate it through cross validation
AllFeatures=[];
AllKinematics=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.Imagined_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.Imagined_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.Imagined_hhmmss{1,j},'Imagined');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
                NeuralFeature=TrialData.NeuralFeatures(:,idx);
                AllFeatures=[AllFeatures, cell2mat(NeuralFeature)];
                AllKinematics=[AllKinematics, TrialData.CursorState(:,idx)];           
        end   
    end
end
% performig regress for individual feature to calculate the C matrices

Day5_Imagiend_Blocks_R2{1,1}='Features';
Day5_Imagiend_Blocks_R2{2,1}='PhaseofDelta';
Day5_Imagiend_Blocks_R2{3,1}='Delta';
Day5_Imagiend_Blocks_R2{4,1}='Theta';
Day5_Imagiend_Blocks_R2{5,1}='Alpha';
Day5_Imagiend_Blocks_R2{6,1}='Beta';
Day5_Imagiend_Blocks_R2{7,1}='LowGamma';
Day5_Imagiend_Blocks_R2{8,1}='HighGamma';

Day5_Imagiend_Blocks_R2{1,2}='Posx';
Day5_Imagiend_Blocks_R2{1,3}='Posy';
Day5_Imagiend_Blocks_R2{1,4}='Velx';
Day5_Imagiend_Blocks_R2{1,5}='Vely';
Day5_Imagiend_Blocks_R2{1,6}='Const';

for Fe=1:7
    Input_dynamics=AllFeatures(128*(Fe-1)+1:(128*(Fe-1)+128),:);
    
    for Kin=1:5
        Output_Dynamics=AllKinematics(Kin,:);
        Figures=0;
        Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
        Day5_Imagiend_Blocks_R2{1+Fe,1+Kin}=Rsq_Pre_output;
        
    end
end


%% Section 3-3: Using CLDA blocks: calculating the new C matrices for individual feature and evaluate it through cross validation
AllFeatures=[];
AllKinematics=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.CLDA_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            
              if length(TrialData.Events)==2 && TrialData.ErrorID==0% it means the cursor start from center and hit the target
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
              
                NeuralFeature=TrialData.NeuralFeatures(:,idx);
                AllFeatures=[AllFeatures, cell2mat(NeuralFeature)];
                AllKinematics=[AllKinematics, TrialData.CursorState(:,idx)];
                
              end       
        end   
    end
end
% performig regress for individual feature to calculate the C matrices

Day5_CLDA_Blocks_R2{1,1}='Features';
Day5_CLDA_Blocks_R2{2,1}='PhaseofDelta';
Day5_CLDA_Blocks_R2{3,1}='Delta';
Day5_CLDA_Blocks_R2{4,1}='Theta';
Day5_CLDA_Blocks_R2{5,1}='Alpha';
Day5_CLDA_Blocks_R2{6,1}='Beta';
Day5_CLDA_Blocks_R2{7,1}='LowGamma';
Day5_CLDA_Blocks_R2{8,1}='HighGamma';

Day5_CLDA_Blocks_R2{1,2}='Posx';
Day5_CLDA_Blocks_R2{1,3}='Posy';
Day5_CLDA_Blocks_R2{1,4}='Velx';
Day5_CLDA_Blocks_R2{1,5}='Vely';
Day5_CLDA_Blocks_R2{1,6}='Const';

for Fe=1:7
    Input_dynamics=AllFeatures(128*(Fe-1)+1:(128*(Fe-1)+128),:);
    
    for Kin=1:5
        Output_Dynamics=AllKinematics(Kin,:);
        Figures=0;
        Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
        Day5_CLDA_Blocks_R2{1+Fe,1+Kin}=Rsq_Pre_output;
        
    end
end

%% Section 3-4: Using Performance blocks: calculating the new C matrices for individual feature and evaluate it through cross validation
AllFeatures=[];
AllKinematics=[];
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.perf_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.perf_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.perf_hhmmss{1,j},'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            
              if length(TrialData.Events)==2 && TrialData.ErrorID==0% it means the cursor start from center and hit the target
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
              
                NeuralFeature=TrialData.NeuralFeatures(:,idx);
                AllFeatures=[AllFeatures, cell2mat(NeuralFeature)];
                AllKinematics=[AllKinematics, TrialData.CursorState(:,idx)];
                
              end       
        end   
    end
end
% performig regress for individual feature to calculate the C matrices

Day5_Fixed_Blocks_R2{1,1}='Features';
Day5_Fixed_Blocks_R2{2,1}='PhaseofDelta';
Day5_Fixed_Blocks_R2{3,1}='Delta';
Day5_Fixed_Blocks_R2{4,1}='Theta';
Day5_Fixed_Blocks_R2{5,1}='Alpha';
Day5_Fixed_Blocks_R2{6,1}='Beta';
Day5_Fixed_Blocks_R2{7,1}='LowGamma';
Day5_Fixed_Blocks_R2{8,1}='HighGamma';

Day5_Fixed_Blocks_R2{1,2}='Posx';
Day5_Fixed_Blocks_R2{1,3}='Posy';
Day5_Fixed_Blocks_R2{1,4}='Velx';
Day5_Fixed_Blocks_R2{1,5}='Vely';
Day5_Fixed_Blocks_R2{1,6}='Const';

for Fe=1:7
    Input_dynamics=AllFeatures(128*(Fe-1)+1:(128*(Fe-1)+128),:);
    
    for Kin=1:5
        Output_Dynamics=AllKinematics(Kin,:);
        Figures=0;
        Rsq_Pre_output=MultipleRegFunc(Input_dynamics',Output_Dynamics',Figures);
        Day5_Fixed_Blocks_R2{1+Fe,1+Kin}=Rsq_Pre_output;
        
    end
end

%% Section 3-5: Make comparison between Imagined results and convergence rate in Fig. 2 

%%%% to make comparison of correlation with convergence rate in Fig 2 paper
%I copied the results of Fig. 2 here:
% STATS: regression of feature decoder angles wrt final decoder
%   reporting regressions of continuous (1st 5 days)
%    feature convergence rates 
% Feature1: int=89.201, slope=-0.067, CI=[-0.072,-0.063]
%   
% Feature2: int=64.369, slope=-0.088, CI=[-0.091,-0.085]
%   
% Feature3: int=68.410, slope=-0.091, CI=[-0.097,-0.085]
%   
% Feature4: int=67.052, slope=-0.022, CI=[-0.032,-0.013]
%   
% Feature5: int=60.445, slope=-0.063, CI=[-0.075,-0.051]
%   
% Feature6: int=59.997, slope=-0.082, CI=[-0.089,-0.076]
%   
% Feature7: int=37.898, slope=-0.082, CI=[-0.087,-0.078]

Avg_Imagined=[...
mean([Day5_Imagiend_Blocks_R2{2, 4},Day5_Imagiend_Blocks_R2{2, 5}]);...
mean([Day5_Imagiend_Blocks_R2{3, 4},Day5_Imagiend_Blocks_R2{3, 5}]);...
mean([Day5_Imagiend_Blocks_R2{4, 4},Day5_Imagiend_Blocks_R2{4, 5}]);...
mean([Day5_Imagiend_Blocks_R2{5, 4},Day5_Imagiend_Blocks_R2{5, 5}]);...
mean([Day5_Imagiend_Blocks_R2{6, 4},Day5_Imagiend_Blocks_R2{6, 5}]);...
mean([Day5_Imagiend_Blocks_R2{7, 4},Day5_Imagiend_Blocks_R2{7, 5}]);...
mean([Day5_Imagiend_Blocks_R2{8, 4},Day5_Imagiend_Blocks_R2{8, 5}])];

Converg_Rates=[-0.067;-0.088;-0.091;-0.022;-0.063;-0.082;-0.082];

[rho,pval]=corr(Avg_Imagined,Converg_Rates);
 

