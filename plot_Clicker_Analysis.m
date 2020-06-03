
% this program is for performing clicker analysis to generate Fig. 5;
% for single clicker A and then stacking theory of multiclicker

clear all
close all
clc 

%% experiment info for single clicker A in GridTask and RadialTask
expts = [];

% Long-Term CLDA in Center-Out using 7 features
% Long-Term CLDA in Center-Out and test performance in square grid 
% CursorControlGridTask
expts(end+1).yymmdd = '20190816';
expts(end).CLDA_hhmmss = {'110820'};
expts(end).perf_GridTask_hhmmss = {'112220','112435','112632','112926','113118','113218',...
    '113506','114851','115108','115346'};

expts(end+1).yymmdd = '20190819';
expts(end).CLDA_hhmmss = {'112813','141607'};
expts(end).perf_GridTask_hhmmss = {'113941','114110','114318','114616','115111','115346',...
    '142226','142626','143335','143651','144113','144505','144816'};

expts(end+1).yymmdd = '20190830';
expts(end).CLDA_hhmmss = {'112710','140606'};
expts(end).perf_GridTask_hhmmss = {'114031','115045','115346','115737','141814','143231',...
    '143449','143810','144326','144640','145048','145735','150124'};

expts(end+1).yymmdd = '20190904';
expts(end).CLDA_hhmmss = {'113405','140826'};
expts(end).perf_GridTask_hhmmss = {'114718','115009','115356','115803','120024',...
    '142452','142649','142912','143715','143850','144117','144534','144902',...
    '145059','145426','145640','145824','150103'};

%%%%%%%%%%%%%%%%%%%%
% please notice from 20190904 to 20191113, we made transition from 
% 7 feature sets to 3 feature sets. We were testing different
% methods(e.g. PCA, Mask) for imprvoing cursor control. Not sure if clicker
% data is good to include here although clicker always uses 6 feature sets.
% Also we made transition from square grid task to radial task

%%%%%%%%%%%%%%%%%%%%%
% Long-Term CLDA in Center-Out using three features: delta, beta, high gamma 
% Testing the performance in Center-Out or Radial tasks, as well

% the following days follow the previous version of Daniel's code
% structures: "Center-Out", "CursorControlGridTask"  
expts(end+1).yymmdd = '20191113';
expts(end).CLDA_hhmmss = {'113926','114736','133000'};
expts(end).perf_CenterOut_hhmmss = {'113926','114736','133000'};
expts(end).perf_RadialTask_hhmmss = {'133705','134023','134920','140150',...
    '140440','141032','141402'};

expts(end+1).yymmdd = '20191115';
expts(end).CLDA_hhmmss = {'111606'};
expts(end).perf_CenterOut_hhmmss = {'111606'};
expts(end).perf_RadialTask_hhmmss = {'144954','145950','150506'};


%% Analysis and stability of single clicker A across days


for day=1:length(expts)
    ClickCounter=0;
    
    expt = expts(day);
    yymmdd = expt.yymmdd;
    
    if day<=4
        BlockNum=length(expt.perf_GridTask_hhmmss);
        Blocks=expt.perf_GridTask_hhmmss;
    elseif day>4
        BlockNum=length(expt.perf_RadialTask_hhmmss);
        Blocks=expt.perf_RadialTask_hhmmss;
    end
    
    for Block=1:BlockNum
        hhmmss=Blocks{1,Block};
        
        datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
            'GangulyServer','CursorControlGridTask',yymmdd,hhmmss,'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        
        for Trial=1:length(files)
            load(fullfile(datadir,files(Trial).name));
            % check the clicker is on and also subject click
            % successfully
            if isfield(TrialData,'ClickerState')&& (TrialData.TargetID==TrialData.SelectedTargetID)
                ClickCounter=ClickCounter+1;
                NeuralData(day).AllClicks(ClickCounter)={cell2mat(TrialData.NeuralFeatures)};
                Threshholds(day).AllClicks(ClickCounter)=TrialData.Params.DecisionBoundary;
                Sample_Freq(day).AllClicks(ClickCounter)=TrialData.Params.UpdateRate;
                
            end
            
        end
        
    end
    
end

%% Generating the clickerspace for all days 
load('clicker_svm_mdl_819-OnlyOK.mat');

for day=1:length(NeuralData)
    
    for Click=1:length(NeuralData(day).AllClicks)
        
        SingleData=cell2mat(NeuralData(day).AllClicks(Click));
        % change the threshhold to zero
        Values=(model.w)*SingleData(129:end,:);%-(-.2);%Threshholds(day).AllClicks(Click);
        ClickerSpace(day).AllClicks(Click)={Values};
        
    end
    
end

% plot an example day to see
for day=1
    for Click=1:length(NeuralData(day).AllClicks)
        SingleData=cell2mat(NeuralData(day).AllClicks(Click));
        Values=(model.w)*SingleData(129:end,:);
        figure;
        stem(Values(1:end-1),'linewidth',1.5)
        hold on
        stem(length(Values)+1,Values(end),'r','linewidth',2)
        hold on
        %from classifer
        hline (Threshholds(day).AllClicks(Click),'m','linewidth',1)
    end
end

% choose one plot for paper

for day=1
    for Click=65
        SingleData=cell2mat(NeuralData(day).AllClicks(Click));
        Values=(model.w)*SingleData(129:end,:)-Threshholds(day).AllClicks(Click);
        
        figure('position',[400 400 500 150]);
        stem(Values(1:end-2),'linewidth',1.5)
        hold on
        stem(length(Values),Values(end-1),'r','linewidth',2)
        hold on
        %from classifer
        hline (0,'m','linewidth',1)
        box off
        ylim([-1 3])
        xlim([0 18])
        yticks([-1,0,1,2,3])
        % sampling rate is 5 so,
        xticks([0 5 10 15])
        xticklabels({'0','1','2','3'})
        set(gca,'FontSize',20)
        
    end
end



%% plotting the clicker space only around click bin. 

for day=1:length(ClickerSpace)
    
    for Click=1:length(ClickerSpace(day).AllClicks)
        
        if Sample_Freq(day).AllClicks(Click)==5
            %upsample the data
            LowRate=ClickerSpace(day).AllClicks{1,Click};
            x=1:length(LowRate);
            xq=0.5:0.5:length(LowRate);
            Value=interp1(x,LowRate,xq);
            
        else
            Value=ClickerSpace(day).AllClicks{1,Click};
            
        end
        
        if size(Value,2)<=21
            Value=[NaN(1,21),Value];
        end
        
        ClickSpaceBin(day).Values(Click,:)=Value(end-21:end);
    end
    
end

for day=1:length(ClickSpaceBin)
    figure('position',[400 200 1000 500]);
    plot(nanmean(ClickSpaceBin(day).Values,1),'b','linewidth',2)
    hold on
    hline (-0.0,'r','linewidth',1)
    HighQualityFigs(['ClickerSpace_Single_Day',num2str(day)])

end


%% Generate the figures with confidence interval across 6 days of clicker A

figure('position',[400 200 600 150]);
for day=1:6
    % plotting with confidence interval
    ValuesA=ClickSpaceBin(day).Values;
    N_A = size(ValuesA,1); % Number of ‘Experiments’ In Data Set
    Mean_A = nanmean(ValuesA); % Mean Of All Experiments At Each Value Of ‘x’
    SEM_A = nanstd(ValuesA)/sqrt(N_A); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], N_A-1);% Calculate 95% Probability Intervals Of t-Distribution
    CI95_A = bsxfun(@times, SEM_A , CI95(:)); %Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    
    % figure
    ciplot(Mean_A+CI95_A(1,:),Mean_A-CI95_A(1,:),30*(day-1)+1:1:(day*30-8),'b')
    %alpha(.4)
    hold on
    hline (-0,'r','linewidth',1)
    set(gca,'FontSize',12)
    %ylabel('Averaged Clicker State Value','fontsize',14)
    %ylim([-1,1])
    %xlabel('Time (mm) relative to click','fontsize',14)
    %xlim([-900 100])
    %xticks([-900 -600 -300 0])
    box off
    %title('Averaged clicker state values aross 4 blocks for Class A')
    %HighQualityFigs('ClickerSpaceA')
    hold on 
    
end 
ylim([-1,2])
yticks([-1 0 1 2])
xlim([1 175])
xticks([1 19, 31 50, 61 81,91 111, 121 141, 151 171])
xticklabels({'-2(s)','0','-2(s)','0','-2(s)','0','-2(s)','0','-2(s)','0','-2(s)','0'})


%% for next section related to two neural clickers 
clear all 
close all
clc
%% experiment info for multiclikcer
expts = [];

expts(end+1).yymmdd = '20191206'; % last day in fig 1
expts(end).hhmmss = {'112558','113508','114445','134919'};


%% Measurement of success rate to select the block for further anlaysis

for day=1:length(expts)
    
    expt = expts(day);
    yymmdd = expt.yymmdd;
    
    for session=1:length(expt.hhmmss)
        hhmmss = expt.hhmmss{session};
        
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer',yymmdd,'RadialTypingMultiClick',hhmmss,'BCI_Fixed');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        
        Fails=0;
        for Trial=1:length(datafiles)
            load(fullfile(datadir,datafiles(Trial).name))
            
            if ~strcmp(TrialData.TargetCharStr,TrialData.SelectedTargetCharStr)%strcmp(TrialData.ErrorStr,'WrongTarget')
                Fails=Fails+1;
                
            end
            
        end
       SuccessRates(session)=(length(datafiles)-Fails)/length(datafiles)*100
         
    end
    
end

% SuccessRates: 45.0000   60.0000   73.3333   58.3333

%% Plotting the projection of neural data to Clicker space for all duration of trials

load('multiclick_svm_kick_vs_ok.mat');

for day=1:length(expts)
    
    expt = expts(day);
    yymmdd = expt.yymmdd;
    
    for session=1:length(expt.hhmmss)
        hhmmss = expt.hhmmss{session};
        
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer',yymmdd,'RadialTypingMultiClick',hhmmss,'BCI_Fixed');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        
        ClickA=0;
        ClickB=0;
        ClickAB=0;
        for Trial=1:length(datafiles)
            load(fullfile(datadir,datafiles(Trial).name))
            
            if strcmp(TrialData.TargetCharStr,TrialData.SelectedTargetCharStr)
                ClickAB=ClickAB+1;
                
                if strcmp(TrialData.TargetCharStr,'A')
                    ClickA=ClickA+1;
                    NeuralData(session).ClassA(ClickA)={cell2mat(TrialData.NeuralFeatures)};
                    Data_A=cell2mat(NeuralData(session).ClassA(ClickA));
                    ClickerSpace(session).ClassA(ClickA)={(model.w)*Data_A(129:end,:)};
                    
                elseif strcmp(TrialData.TargetCharStr,'B')
                    ClickB=ClickB+1;
                    NeuralData(session).ClassB(ClickB)={cell2mat(TrialData.NeuralFeatures)};
                    Data_B=cell2mat(NeuralData(session).ClassB(ClickB));
                    ClickerSpace(session).ClassB(ClickB)={(model.w)*Data_B(129:end,:)};
                    
                    
                end
                
                VelMagAll(session).VelMag(ClickAB)={sqrt(TrialData.CursorState(3,:).^2+TrialData.CursorState(4,:).^2)};
                
            end
            
        end
       
         
    end
    
end

for session=1:length(ClickerSpace)
    
    PlotNum=max([length(ClickerSpace(session).ClassA),length(ClickerSpace(session).ClassB)]);
    
    figure('position',[400 200 1500 1100]);
    suptitle(['Session:',expt.hhmmss{session},'; Clicker space values for correct selections: Class A (left column); Class B (right column)'])
    
    for i=1:length(ClickerSpace(session).ClassA)
        subplot(PlotNum,2,i*2-1)
        ValuesA=ClickerSpace(session).ClassA{1,i};
        stem(ValuesA(1:end-1),'linewidth',1.5)
        hold on
        stem(length(ValuesA)+1,ValuesA(end),'r','linewidth',2)
        hold on
        %from classifer
        hline (-0.4,'m','linewidth',1)
        
    end
    
    for i=1:length(ClickerSpace(session).ClassB)
        
        subplot(PlotNum,2,i*2)
        ValuesB=ClickerSpace(session).ClassB{1,i};
        stem(ValuesB(1:end-1),'linewidth',1.5)
        hold on
        stem(length(ValuesB)+1,ValuesB(end),'r','linewidth',2)
        hold on
        %from classifer
        hline (0.2,'m','linewidth',1)
        
    end

%HighQualityFigs(['ClickerSpace_Session',expt.hhmmss{session}])    
end


%% plotting the clicker space only around click bin. 

CounterA=0;
CounterB=0;
ValuesA=[];
ValuesB=[];
for session=1:length(ClickerSpace)
    
    for i=1:length(ClickerSpace(session).ClassA)
        CounterA=CounterA+1;
        Value=ClickerSpace(session).ClassA{1,i};
        if size(Value,2)<=21
            Value=[NaN(1,21),Value];
        end
        ValuesA(CounterA,:)=Value(end-21:end);   
    end
    
    for i=1:length(ClickerSpace(session).ClassB)
        CounterB=CounterB+1;
        Value=ClickerSpace(session).ClassB{1,i};
        if size(Value,2)<=21
            Value=[NaN(1,21),Value];
        end
        ValuesB(CounterB,:)=Value(end-21:end);   
    end

end

figure('position',[400 200 1000 500]);
plot(nanmean(ValuesA,1),'b','linewidth',2)
hold on
hline (-0.4,'g','linewidth',1)
hold on
plot(nanmean(ValuesB,1),'r','linewidth',2)
hold on 
hline (0.2,'m','linewidth',1)
legend('Class A','Boundary A','Class B','Boundary B','Location','Northwest')
set(gca,'FontSize',18)

ylabel('Averaged Clicker State Value','fontsize',14)
ylim([-1,1])
xlabel('Bins before Click','fontsize',14)
xlim([1 22])
title('Averaged clicker state values aross 4 blocks for Class A and Class B')
%HighQualityFigs('ClickerSpaceAB')


% plotting with confidence interval
N_A = size(ValuesA,1); % Number of ‘Experiments’ In Data Set
Mean_A = nanmean(ValuesA); % Mean Of All Experiments At Each Value Of ‘x’
SEM_A = nanstd(ValuesA)/sqrt(N_A); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N_A-1);% Calculate 95% Probability Intervals Of t-Distribution
CI95_A = bsxfun(@times, SEM_A , CI95(:)); %Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

N_B = size(ValuesB,1); % Number of ‘Experiments’ In Data Set
Mean_B = nanmean(ValuesB); % Mean Of All Experiments At Each Value Of ‘x’
SEM_B = nanstd(ValuesB)/sqrt(N_B); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N_B-1);% Calculate 95% Probability Intervals Of t-Distribution
CI95_B = bsxfun(@times, SEM_B , CI95(:)); %Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

% figure
figure('position',[400 200 400 150]);
ciplot(Mean_A+CI95_A(1,:),Mean_A-CI95_A(1,:),1:1:22,'b')
%alpha(.4)
hold on
hline (-0.52,'r','linewidth',1)
set(gca,'FontSize',18)
%ylabel('Averaged Clicker State Value','fontsize',14)
ylim([-1,1])
%xlabel('Time (mm) relative to click','fontsize',14)
xlim([1 22])
xticks([2,12,22])
xticklabels({'-2','-1','0'})
box off
%title('Averaged clicker state values aross 4 blocks for Class A')
%HighQualityFigs('ClickerSpaceA')

figure('position',[400 200 400 150]);
ciplot(Mean_B+CI95_B(1,:),Mean_B-CI95_B(1,:),1:1:22,'b')
%alpha(.4)
hold on
hline (0.2,'r','linewidth',1)
set(gca,'FontSize',18)
%ylabel('Averaged Clicker State Value','fontsize',14)
ylim([-1,1])
%xlabel('Time (mm) relative to click','fontsize',14)
xlim([1 22])
xticks([2,12,22])
xticklabels({'-2','-1','0'})
box off
%title('Averaged clicker state values aross 4 blocks for Class B')
%HighQualityFigs('ClickerSpaceB')


%% Investigate the neural features values at the moment of click:  Method 1
CounterA=0;
CounterB=0;
FeatureA=[];
FeatureB=[];
% for all features
for session=1:length(NeuralData)
    
    for i=1:length(NeuralData(session).ClassA)
        CounterA=CounterA+1;
        FeatureA(:,CounterA)=NeuralData(session).ClassA{1,i}(128+1:end,end);
    end
    
    for i=1:length(NeuralData(session).ClassB)
        CounterB=CounterB+1;
        FeatureB(:,CounterB)=NeuralData(session).ClassB{1,i}(128+1:end,end);
    end
    
end

ch_layout = [
    96	84	76	95	70	82	77	87	74	93	66	89	86	94	91	79
    92	65	85	83	68	75	78	81	72	69	88	71	80	73	90	67
    62	37	56	48	43	44	60	33	49	64	58	59	63	61	51	34
    45	53	55	52	35	57	38	50	54	39	47	42	36	40	46	41
    19	2	10	21	30	23	17	28	18	1	8	15	32	27	9	3
    24	13	6	4	7	16	22	5	20	14	11	12	29	26	31	25
    124	126	128	119	110	113	111	122	117	125	112	98	104	116	103	106
    102	109	99	101	121	127	105	120	107	123	118	114	108	115	100	97];

feature_str = {'|\delta|','|\theta|','|\alpha|',...
    '|\beta|','|\gamma_1|','|\gamma_2|'};

MeanA_map=mean(FeatureA,2);
MeanB_map=mean(FeatureB,2);
FeatureA_map=[];
FeatureB_map=[];
Diff_AB_map=[];

figure('Position',[673 100 540 1152]);
for f=1:6
    for ch=1:128
        [r,c] = find(ch_layout==ch);
        FeatureA_map(r,c) = (MeanA_map(ch+128*(f-1),1));
        FeatureB_map(r,c) = (MeanB_map(ch+128*(f-1),1));
        Diff_AB_map(r,c)=(MeanA_map(ch+128*(f-1),1)-MeanB_map(ch+128*(f-1),1));
        
    end
    
    subplot(6,3,3*f-2);
    imagesc(FeatureA_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},' feature: Clicker A'])
    
    subplot(6,3,3*f-1);
    imagesc(FeatureB_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},' feature: Clicker B'])
    
    subplot(6,3,3*f);
    imagesc(Diff_AB_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},' feature: Clicker A-B'])
    
    
end

%HighQualityFigs('ClickerMapAB_1')



%% Investigate the neural features values*SVM Weight at the moment of click:  Method 2

CounterA=0;
CounterB=0;
FeatureA=[];
FeatureB=[];
% for high gamma only
for session=1:length(NeuralData)
    
    for i=1:length(NeuralData(session).ClassA)
        CounterA=CounterA+1;
        FeatureA(:,CounterA)=NeuralData(session).ClassA{1,i}(128+1:end,end);
    end
    
    for i=1:length(NeuralData(session).ClassB)
        CounterB=CounterB+1;
        FeatureB(:,CounterB)=NeuralData(session).ClassB{1,i}(128+1:end,end);
    end
    
end

ch_layout = [
    96	84	76	95	70	82	77	87	74	93	66	89	86	94	91	79
    92	65	85	83	68	75	78	81	72	69	88	71	80	73	90	67
    62	37	56	48	43	44	60	33	49	64	58	59	63	61	51	34
    45	53	55	52	35	57	38	50	54	39	47	42	36	40	46	41
    19	2	10	21	30	23	17	28	18	1	8	15	32	27	9	3
    24	13	6	4	7	16	22	5	20	14	11	12	29	26	31	25
    124	126	128	119	110	113	111	122	117	125	112	98	104	116	103	106
    102	109	99	101	121	127	105	120	107	123	118	114	108	115	100	97];

feature_str = {'|\delta|','|\theta|','|\alpha|',...
    '|\beta|','|\gamma_1|','|\gamma_2|'};

MeanA_map=mean(FeatureA,2);
MeanB_map=mean(FeatureB,2);
FeatureA_map=[];
FeatureB_map=[];
Diff_AB_map=[];

load('multiclick_svm_kick_vs_ok.mat');

figure('Position',[673 100 540 1152]);
for f=1:6
    for ch=1:128
        [r,c] = find(ch_layout==ch);
        FeatureA_map(r,c) = model.w(ch+128*(f-1))*(MeanA_map(ch+128*(f-1),1));
        FeatureB_map(r,c) = model.w(ch+128*(f-1))*(MeanB_map(ch+128*(f-1),1));
        Diff_AB_map(r,c)=model.w(ch+128*(f-1))*(MeanA_map(ch+128*(f-1),1)-MeanB_map(ch+128*(f-1),1));
        
    end
    
    subplot(6,3,3*f-2);
    imagesc(FeatureA_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},'Fe*W(SVM): Clicker A'])
    
    subplot(6,3,3*f-1);
    imagesc(FeatureB_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},'Fe*W(SVM): Clicker B'])
    
    subplot(6,3,3*f);
    imagesc(Diff_AB_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    %colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},' The diff: Clicker A-B'])
    
    
end

%HighQualityFigs('ClickerMapAB_2')

%% Investigate the weights of SVM at the moment of click and the diff with clicker A:  Method 3

clear model
load('clicker_svm_mdl_819-OnlyOK.mat');
W_A=model.w;
load('multiclick_svm_kick_vs_ok.mat');
W_B=model.w;

% figure('Position',[673 50 10 950]);
% for f=1:6
%     SVM_map=[];
%     for ch=1:128
%         [r,c] = find(ch_layout==ch);
%         SVM_map(r,c) = abs(model.w(ch+128*(f-1)));
%     end
%     
%     subplot(6,1,f);
%     imagesc(SVM_map);
%     colorbar('southoutside')
%     set(gca,'xtick',[],'ytick',[])
%     colormap(brewermap(100,'YlOrRd'))
%     title([feature_str{f},'; W(SVM)'])
%     
% end

%HighQualityFigs('ClickerMapSVM_AB')

% only for clicker B in high gamma: for figure inside paper
figure('Position',[673 200 200 200]);
for f=6
    SVM_map=[];
    for ch=1:128
        [r,c] = find(ch_layout==ch);
        SVM_map(r,c) = abs(W_B(ch+128*(f-1)));
    end
    
    %subplot(6,1,f);
    imagesc(SVM_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},'; WB(SVM)'])
    
end


% the diff with Clicker A in high gamma: for figure inside paper
figure('Position',[673 200 200 200]);
for f=6
    SVM_map=[];
    for ch=1:128
        [r,c] = find(ch_layout==ch);
        SVM_map(r,c) = abs(abs(W_B(ch+128*(f-1)))-abs(W_A(ch+128*(f-1))));
    end
    
    %subplot(6,1,f);
    imagesc(SVM_map);
    colorbar('southoutside')
    set(gca,'xtick',[],'ytick',[])
    colormap(brewermap(100,'YlOrRd'))
    title([feature_str{f},'; WB-WA(SVM)'])
    
end



%% calculated the speed of cursor during movement versus during clicking

Speed_Cursor=[];
Speed_Click=[];

for session=1:length(VelMagAll)
    
    for i=1:length(VelMagAll(session).VelMag)
        Speed_Cursor=[Speed_Cursor, VelMagAll(session).VelMag{1,i}(1,1:end-1)];
        Speed_Click=[Speed_Click,VelMagAll(session).VelMag{1,i}(1,end)];
    end
    
    
end

figure('Position',[673 100 360 360]);
hist(Speed_Cursor)
xlabel('Speed (Scaled)')
xlim([0 40])
xticks([10 20 30 40])
ylabel('Number of Bins')
ylim([0,800])
yticks([0 200 400 600 800])
title('Dist. of cursor speed: During Travelling')
set(gca,'FontSize',18)
box off
HighQualityFigsSVG('R4_1')

figure('Position',[673 100 360 360]);
hist(Speed_Click)
xlabel('Speed (Scaled)')
xlim([0 40])
xticks([10 20 30 40])
ylabel('Number of Bins')
ylim([0,30])
yticks([0 10 20 30])
title('Dist. of cursor speed: During click')
set(gca,'FontSize',18)
box off
HighQualityFigsSVG('R4_2')

%HighQualityFigs('ClickerSpeed_AB')

 