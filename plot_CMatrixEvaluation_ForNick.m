% This is the analysis for reviewers comments regarding to individual
% feature effect even before ltCLDA or on first day of ltCLDA:

%% Section 1: using the day 
% features

clear all,
close all,
clc

%% Sestion 1-1: laoding the first day of ltCLDA
% experiment info:
% Long-Term CLDA in Center-Out and test performance in Center-Out 
expts = [];
expts(end+1).yymmdd = '20190521';
expts(end).Imagined_hhmmss = {'133731'};%{'110055','113915','114657'};
expts(end).CLDA_hhmmss = {'135008','141438'};
expts(end).perf_hhmmss = {'135943','141438'};

%% Sestion 1-1: laoding the last day of imagined blocks

% experiment info:
% CLDA in Center-Out and test performance in Center-Out 
expts = [];
expts(end+1).yymmdd = '20190515';
expts(end).Imagined_hhmmss = {'132722'};
expts(end).CLDA_hhmmss = {'110735'};
expts(end).perf_hhmmss = {'112302','113447'};

%% Section 2-2: calculating the new C matrices for all feature and evaluate it through cross validation
% also adding imagined blocks/days to calculation to see the if C got better

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
        
       
        %for performig regress for all features together per session to
        %calculate the C matrix
        Input_dynamics=PerSession_Features;
        for Kin=1:5
            Output_Dynamics=PerSession_Kinematics(Kin,:);
            Figures=0;
            CMatrix=MultipleRegFunc_Nick(Input_dynamics',Output_Dynamics',Figures);
            % the row is feature the column is the kin
            CMatrix_all(i).PerSession(j).AllFe(Kin)={CMatrix};
                       
        end
        
    end
      
end

Kin_Parameter={'Posx','Posy','Velx','Vely','Constant'};

save('CMatrix_Day0.mat','CMatrix_all')

