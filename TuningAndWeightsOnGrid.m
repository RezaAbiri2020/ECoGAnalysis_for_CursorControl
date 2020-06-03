

%% Plotting tuning and weights on the grid for fixed blocks (For both center reset on/off/discrete/continuous control sessions).
% the program cab be run for imagined or adapt blocks as well
clear all
close all
NCh=128;

% Potential days for analysis:
% 20190514; day1
% 20190515; day2
% 20190521; day3
% 20190524; day4
% 20190529; day5
% 20190531; day6

%% including the selected Fixed sessions for each day for analysis through
if 0
    % for day1: April26
    day1='../Data\20190426\GangulyServer\Center-Out\20190426\';
    Bravo1(1).Days={[day1,'115111\BCI_Fixed'],...
        [day1,'144334\BCI_Fixed']};
    
    % for day2: April29
    day2='../Data\20190429\GangulyServer\Center-Out\20190429\';
    Bravo1(2).Days={[day2,'130907\BCI_Fixed'],...
        [day2,'134613\BCI_Fixed'],...
        [day2,'135428\BCI_Fixed'],...
        [day2,'140100\BCI_Fixed']};
    
    % for day3: May01
    day3='../Data\20190501\GangulyServer\Center-Out\20190501\';
    Bravo1(3).Days={[day3,'105345\BCI_Fixed'],...
        [day3,'112049\BCI_Fixed'],...
        [day3,'135512\BCI_Fixed']};
end

if 0
    % for day1: May14
    day1='../Data\20190514\GangulyServer\Center-Out\20190514\';
    Bravo1(1).Days={[day1,'135927\BCI_Fixed'],...
        [day1,'140808\BCI_Fixed'],...
        [day1,'141731\BCI_Fixed'],...
        [day1,'142704\BCI_Fixed']};
    
    % for day2: May15
    day2='../Data\20190515\GangulyServer\Center-Out\20190515\';
    Bravo1(2).Days={[day2,'112302\BCI_Fixed'],...
        [day2,'113447\BCI_Fixed'],...
        [day2,'135831\BCI_Fixed'],...
        [day2,'141134\BCI_Fixed']...
        [day2,'141859\BCI_Fixed']};
    
    % for day3: May21
    day3='../Data\20190521\GangulyServer\20190521\';
    Bravo1(3).Days={[day3,'135943\BCI_Fixed'],...
        [day3,'141438\BCI_Fixed']};
end

if 1
    % for day4: May24
    day4='../Data\20190524\GangulyServer\Center-Out\20190524\';
    Bravo1(4).Days={[day4,'111731\BCI_Fixed'],...
        [day4,'112353\BCI_Fixed'],...
        [day4,'134653\BCI_Fixed'],...
        [day4,'135957\BCI_Fixed']};
    
    % for day5: May29
    day5='../Data\20190529\GangulyServer\Center-Out\20190529\';
    Bravo1(5).Days={[day5,'111528\BCI_Fixed'],...
        [day5,'114050\BCI_Fixed'],...
        [day5,'135644\BCI_Fixed'],...
        [day5,'141529\BCI_Fixed']};
    
    % for day6: May31
    day6='../Data\20190531\GangulyServer\Center-Out\20190531\';
    Bravo1(6).Days={[day6,'105048\BCI_Fixed'],...
        [day6,'110703\BCI_Fixed'],...
        [day6,'112444\BCI_Fixed'],...
        [day6,'133517\BCI_Fixed'],...
        [day6,'140204\BCI_Fixed'],...
        [day6,'141319\BCI_Fixed']};
end

%% including the selected adapt sessions for each day for analysis through
if 1
    % for day1: May14
    day1='../Data\20190514\GangulyServer\Center-Out\20190514\';
    Bravo1(1).Days={[day1,'112836\BCI_CLDA'],...
        [day1,'133050\BCI_CLDA']};
    
    % for day2: May15
    day2='../Data\20190515\GangulyServer\Center-Out\20190515\';
    Bravo1(2).Days={[day2,'110735\BCI_CLDA'],...
        [day2,'134017\BCI_CLDA']};
    
    % for day3: May21
    day3='../Data\20190521\GangulyServer\20190521\';
    Bravo1(3).Days={[day3,'135008\BCI_CLDA'],...
        [day3,'141438\BCI_CLDA']};
end

if 1
    % for day4: May24
    day4='../Data\20190524\GangulyServer\Center-Out\20190524\';
    Bravo1(4).Days={[day4,'110219\BCI_CLDA'],...
        [day4,'133313\BCI_CLDA']};
    
    % for day5: May29
    day5='../Data\20190529\GangulyServer\Center-Out\20190529\';
    Bravo1(5).Days={[day5,'105609\BCI_CLDA'],...
        [day5,'113247\BCI_CLDA'],...
        [day5,'133450\BCI_CLDA']};
    
    % for day6: May31
    day6='../Data\20190531\GangulyServer\Center-Out\20190531\';
    Bravo1(6).Days={[day6,'102944\BCI_CLDA'],...
        [day6,'111946\BCI_CLDA'],...
        [day6,'132244\BCI_CLDA'],...
        [day6,'135046\BCI_CLDA']};
end



%% including the imagined sessions for each day for analysis through

if 1
    % for day1: May14
    day1='../Data\20190514\GangulyServer\Center-Out\20190514\';
    Bravo1(1).Days={[day1,'110200\Imagined'],...
        [day1,'111822\Imagined']};
    
    % for day2: May15
    day2='../Data\20190515\GangulyServer\Center-Out\20190515\';
    Bravo1(2).Days={[day2,'105804\Imagined'],...
        [day2,'132722\Imagined']};
    
    % for day3: May21
    day3='../Data\20190521\GangulyServer\20190521\';
    Bravo1(3).Days={[day3,'110055\Imagined'],...
        [day3,'113915\Imagined'],...
        [day3,'114657\Imagined'],...
        [day3,'133731\Imagined']};
end



%% Extracting the necessary data for analysis

for day=1:length(Bravo1)
    
    for session=1:length(Bravo1(day).Days)
        
        % Loading the related data set...
        datadir=char(Bravo1(day).Days(session));
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        
        % extracting the weight values for vx and vy
        % load one trial to get the weights out of it
        load(fullfile(datadir,datafiles(end).name))
        
        %number of features
        NFeatures=TrialData.Params.NumFeatures;
        %Check if zscore happened for all features
        zfeatures=TrialData.Params.ZscoreFeaturesFlag;
        
        Coeff_Features=TrialData.KalmanFilter{1,end}.C;
        
        % third column for vx and forth column for vy
        %plot(Coeff_Features(:,3))
        %plot(Coeff_Features(:,4))
        
        % for vx
        Vx_Coeffs=reshape(Coeff_Features(:,3),NCh,NFeatures);
        % for vy
        Vy_Coeffs=reshape(Coeff_Features(:,4),NCh,NFeatures);
        
        % organizing the weight for plotting
        ch_layout = [
            91	84	67	90	70	79	88	69	92	83	65	89	87	86	94	82
            66	93	78	95	76	75	85	73	68	80	74	72	96	71	77	81
            60	37	42	50	56	54	49	40	43	35	45	63	47	46	58	55
            53	57	33	48	39	51	41	34	64	52	62	38	36	44	61	59
            8	26	29	28	9	5	13	20	11	23	16	22	27	4	3	31
            7	21	15	24	25	1	2	32	14	12	30	19	18	17	6	10
            110	125	111	115	103	117	100	123	113	119	118	98	101	105	116	99
            107	112	97	128	121	124	108	109	127	126	106	122	114	120	104	102];
        ch_layout=ch_layout';
        ch_layout=ch_layout(:);
        Vx_Coeffs_Brain=zeros(NCh,NFeatures);
        Vy_Coeffs_Brain=zeros(NCh,NFeatures);
        
        for j=1:NCh
            Vx_Coeffs_Brain(j,:)=Vx_Coeffs(ch_layout(j),:);
            Vy_Coeffs_Brain(j,:)=Vy_Coeffs(ch_layout(j),:);
        end
        
        % I have Coeffs_Brain for this structure for each feature:
        % 1 2 3 ..
        % 17 18....
        % ........128
        
        AllWeights(day).Days(session).Sessions.VxPerFeature=Vx_Coeffs_Brain;
        AllWeights(day).Days(session).Sessions.VyPerFeature=Vy_Coeffs_Brain;
        
        % Extracting the variables values for all targets
        
        % TargetID: 1: rightTarget (0 degree), then rotating clockwise for other targets ID
        % trial counter per target
        trial_counter1=0;
        trial_counter2=0;
        trial_counter3=0;
        trial_counter4=0;
        trial_counter5=0;
        trial_counter6=0;
        trial_counter7=0;
        trial_counter8=0;
        
        for trial=[1:length(datafiles)]
            %trial=1:NBlocks*8
            %1:length(datafiles)
            %(length(datafiles)-(NBlocks*8-1)):length(datafiles)
            % load data and classifying based on TargetID
            load(fullfile(datadir,datafiles(trial).name))
            
            if length(TrialData.Events)==2
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
                
                % Finding the targetID
                TargetID=TrialData.TargetID;
                % Finding the duration of neurofeedback
                
                % For different blocks; the extracted neural feature data will be different
                if contains(datadir,'Imagined')
                    Nbin=length(idx);
                    
                elseif contains(datadir,'BCI_Fixed') || contains(datadir,'BCI_CLDA') % one feature value per bin
                    Nbin_Hit=length(idx);
                    TimeForAnalysis=2; % second
                    
                    if TrialData.NeuralSamps(5)<150 % binsize is 100ms
                        Nbin=(10*TimeForAnalysis);
                        
                    elseif TrialData.NeuralSamps(5)>150 % binsize is 200ms
                        Nbin=(5*TimeForAnalysis);
                    end
                    
                end
                
                % Extracting data: each row will be for a target with corresponding
                % trials inside that
                if TargetID==1
                    trial_counter1=trial_counter1+1;
                    CursorStates(TargetID).Targets(trial_counter1).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter1).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter1).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter1).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter1).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==2
                    trial_counter2=trial_counter2+1;
                    CursorStates(TargetID).Targets(trial_counter2).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter2).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter2).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter2).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter2).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==3
                    trial_counter3=trial_counter3+1;
                    CursorStates(TargetID).Targets(trial_counter3).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter3).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter3).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter3).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter3).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==4
                    trial_counter4=trial_counter4+1;
                    CursorStates(TargetID).Targets(trial_counter4).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter4).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter4).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter4).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter4).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==5
                    trial_counter5=trial_counter5+1;
                    CursorStates(TargetID).Targets(trial_counter5).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter5).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter5).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter5).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter5).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==6
                    trial_counter6=trial_counter6+1;
                    CursorStates(TargetID).Targets(trial_counter6).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter6).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter6).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter6).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter6).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==7
                    trial_counter7=trial_counter7+1;
                    CursorStates(TargetID).Targets(trial_counter7).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter7).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter7).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter7).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter7).Trials=TrialData.NeuralSamps;
                    
                elseif TargetID==8
                    trial_counter8=trial_counter8+1;
                    CursorStates(TargetID).Targets(trial_counter8).Trials=TrialData.CursorState(:,idx);
                    ECoGs=TrialData.BroadbandData(:,idx);
                    ECoG(TargetID).Targets(trial_counter8).Trials=cell2mat(ECoGs');
                    SampleTimes(TargetID).Targets(trial_counter8).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter8).Trials=cell2mat(NeuralFeature);
                    SampleNumPoints(TargetID).Targets(trial_counter8).Trials=TrialData.NeuralSamps;
                    
                end
                
            end
        end
        
        % Fining the values for Tuning ..
        
        %Feature 1: Phase of delta
        %Feature 2: delta power
        %Feature 3: theta power
        %Feature 4: alpha power
        %Feature 5: beta  power
        %Feature 6: gamma1 power
        %Feature 7: gamma2 power
        
        NFeatures=7;
        FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};
        
        % calculating all features, their means, STD, for all targets and save in a matrix
        for TargetID=1:length(NeuralFeatures) %number of targets
            
            for trial=1:length(NeuralFeatures(TargetID).Targets)
                Features_PerTgTr=nanmean(NeuralFeatures(TargetID).Targets(trial).Trials,2);
                
                ReshapeFeatures=reshape(Features_PerTgTr,NCh,NFeatures);
                
                ch_layout = [
                    91	84	67	90	70	79	88	69	92	83	65	89	87	86	94	82
                    66	93	78	95	76	75	85	73	68	80	74	72	96	71	77	81
                    60	37	42	50	56	54	49	40	43	35	45	63	47	46	58	55
                    53	57	33	48	39	51	41	34	64	52	62	38	36	44	61	59
                    8	26	29	28	9	5	13	20	11	23	16	22	27	4	3	31
                    7	21	15	24	25	1	2	32	14	12	30	19	18	17	6	10
                    110	125	111	115	103	117	100	123	113	119	118	98	101	105	116	99
                    107	112	97	128	121	124	108	109	127	126	106	122	114	120	104	102];
                
                ch_layout=ch_layout';
                ch_layout=ch_layout(:);
                Features_Brain=zeros(NCh,NFeatures);
                
                for j=1:NCh
                    Features_Brain(j,:)=ReshapeFeatures(ch_layout(j),:);
                end
                
                % I have Features_Brain for this structure for each feature:
                % 1 2 3 ..
                % 17 18....
                % ........128
                
                % generati one value per trial per target
                All_Features(TargetID).Targets(trial).Trials=Features_Brain;
                %All_Features(TargetID).Targets(trial).Trials=abs(Features_Brain);
                
            end
            
        end
        
        % calculating the regress/tuning for all targets
        % f=B*X
        % by assuming the chosen tested targets are symmetric in workspace
        section=360/length(NeuralFeatures);
        angles=[0:-section:-(360-section)]; % corresponding to each trial_counter
        
        All_Counters=[trial_counter1,trial_counter2,trial_counter3,trial_counter4,trial_counter5,trial_counter6...
            ,trial_counter7,trial_counter8];
        X=[];
        for count=1:length(All_Counters)
            if All_Counters(count)~=0
                X1=[ones(1,All_Counters(count));ones(1,All_Counters(count))*sind(angles(count));ones(1,All_Counters(count))*cosd(angles(count))];
                X=[X,X1];
            end
            
        end
        
        % finding the values related to tuning
        for Fe=1:NFeatures
            
            for ch=1:NCh
                
                f=[];
                for TargetID=1:length(NeuralFeatures)
                    
                    for trial=1:length(All_Features(TargetID).Targets)
                        f=[f; All_Features(TargetID).Targets(trial).Trials(ch,Fe)];
                    end
                    
                end
                
                [b,bint,r,rint,stats] = regress(f,X');
                PD(ch)=atan2(b(2), b(3)); % radian; to change to degree: /pi*180
                DepthModul(ch)=sqrt(b(2)^2+b(3)^2)/abs(b(1));
                R2(ch)=stats(1);
                
            end
            AllTunings(day).Days(session).Sessions(Fe).Features.PD=PD;
            AllTunings(day).Days(session).Sessions(Fe).Features.DepthModul=DepthModul;
            AllTunings(day).Days(session).Sessions(Fe).Features.R2=R2;
            
            
        end
        
        clear Vx_Coeffs_Brain Vy_Coeffs_Brain;
        clear NeuralFeatures All_Features
        
    end

    
end

%% plotting the tuning 

% The tuning plots should be only for three days for analysis
StartDay=3; %First day of analysis
EndDay=StartDay+2;

for Fe=2 %1:NFeatures
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{Fe},'; Distributions of tunings (FixedBlocks) for each Ch for May24 (Yellow), May29 (green), May31(magenta)'])
    
    for ch=1:NCh
        
        R2_Ch_All=[];
        for day=StartDay:EndDay %length(Bravo1)
            
            PD_Ch_All=[];
            for session=1:length(Bravo1(day).Days)
                
                PD_Ch=AllTunings(day).Days(session).Sessions(Fe).Features.PD(ch);
                PD_Ch_All=[PD_Ch_All;PD_Ch];
                %AllTunings(day).Days(session).Sessions(Fe).Features.DepthModul=DepthModul;
                R2_Ch=AllTunings(day).Days(session).Sessions(Fe).Features.R2(ch);
                R2_Ch_All=[R2_Ch_All;R2_Ch];
                
            end
            
            PD_Ch_mean=mean(PD_Ch_All);
            PD_Ch_std=std(PD_Ch_All);
            level1=PD_Ch_mean-PD_Ch_std;
            level2=PD_Ch_mean+PD_Ch_std;        
            
            subplot(8,16,ch)
            hold on
            alpha = linspace(min(level1,level2),max(level1,level2));
            if day==StartDay
                patch([0 cos(alpha) 0], [0 sin(alpha) 0], 'y','FaceAlpha',.2,'EdgeAlpha',.3)
            elseif day==(StartDay+1)
                patch([0 cos(alpha) 0], [0 sin(alpha) 0], 'g','FaceAlpha',.2,'EdgeAlpha',.3)
            elseif day==(StartDay+2)
                patch([0 cos(alpha) 0], [0 sin(alpha) 0], 'm','FaceAlpha',.2,'EdgeAlpha',.3)
            end
            limit=1;
            ylim([-limit +limit])
            xlim([-limit +limit])
            hold on
            ezplot('x^2+y^2=1');
            axis off
             
        end
        R2_Ch_mean=mean(R2_Ch_All);
        title(['Ch',num2str(ch_layout(ch)),';R2:',num2str(R2_Ch_mean,'%1.2f')]);
       
        
    end
    
    HighQualityFigs([FeatureName{Fe},'_TuningsFixed_Dist3Days(May24-29-31)'])
    
    
end

%% Plotting the weights

for Fe=1:NFeatures
    
    Vx_Fe_Days=[];
    Vy_Fe_Days=[];
    for day=1:length(Bravo1)
        
        Vx_Fe_All=[];
        Vy_Fe_All=[];
        for session=1:length(Bravo1(day).Days)
            
            Vx_Fe=AllWeights(day).Days(session).Sessions.VxPerFeature(:,Fe);
            Vx_Fe_All=[Vx_Fe_All,Vx_Fe];
            
            Vy_Fe=AllWeights(day).Days(session).Sessions.VyPerFeature(:,Fe);
            Vy_Fe_All=[Vy_Fe_All,Vy_Fe];
        end
        
        Vx_Fe_Days(:,day)=mean(Vx_Fe_All,2);
        Vy_Fe_Days(:,day)=mean(Vy_Fe_All,2);
    end
    
    for day=1:length(Bravo1)
        Wx1=reshape(Vx_Fe_Days(:,day),16,8);
        Wy1=reshape(Vy_Fe_Days(:,day),16,8);
        
        aclim=-0.0005;
        bclim=0.0005;
        clims = [aclim bclim];
        
        figure(100);
        Wx2=Wx1';
        imagesc(Wx2,clims);
        colorbar;
        close Figure 100
        % Now make an RGB image that matches display from IMAGESC:
        figure(100);
        Cx = colormap;  % Get the figure's colormap.
        close Figure 100
        Lx = size(Cx,1);
        % Scale the matrix to the range of the map.
        Gsx = round(interp1(linspace(min(Wx2(:)),max(Wx2(:)),Lx),1:Lx,Wx2));
        Hx = reshape(Cx(Gsx,:),[size(Gsx) 3]); % Make RGB image from scaled.
        ColorWeights(day).Vx=Hx;
        
        
        figure(100);
        Wy2=Wy1';
        imagesc(Wy2,clims);
        colorbar;
        close Figure 100
        % Now make an RGB image that matches display from IMAGESC:
        figure(100);
        Cy = colormap;  % Get the figure's colormap.
        close Figure 100
        Ly = size(Cy,1);
        % Scale the matrix to the range of the map.
        Gsy = round(interp1(linspace(min(Wy2(:)),max(Wy2(:)),Ly),1:Ly,Wy2));
        Hy = reshape(Cy(Gsy,:),[size(Gsy) 3]); % Make RGB image from scaled.
        ColorWeights(day).Vy=Hy;
        
    end
    
    
    for xory=1:2
        if xory==1
            figure(Fe),
            set(gcf, 'Position', [100, 100, 2400, 1200]);
            suptitle([FeatureName{Fe},'; Variation of Vx weights (FixedBlocks) in each Ch from left to right corresponding to May14-15-21-24-29-31'])
            
        elseif xory==2
            figure(Fe+NFeatures),
            set(gcf, 'Position', [100, 100, 2400, 1200]);
            suptitle([FeatureName{Fe},'; Variation of Vy weights (FixedBlocks) in each Ch from left to right corresponding to May14-15-21-24-29-31'])
        end
        
        for ch=1:NCh
            row=floor(ch/16)+1;
            if mod(ch,16)~=0
                col=mod(ch,16);
            else
                col=16;
                row=floor(ch/16);
            end
            
            R2_Ch_All=[];
            for day=1:length(Bravo1)
                for session=1:length(Bravo1(day).Days)
                    R2_Ch=AllTunings(day).Days(session).Sessions(Fe).Features.R2(ch);
                    R2_Ch_All=[R2_Ch_All;R2_Ch];
                end
            end
            R2_Ch_mean=mean(R2_Ch_All);
            
            subplot(8,16,ch)
            heightbar=2;
            for day=1:length(Bravo1)
                if xory==1
                    Colors=ColorWeights(day).Vx;
                elseif xory==2
                    Colors=ColorWeights(day).Vy;
                end
                bar(day,heightbar,'FaceColor',Colors(row,col,:))
                axis off
                hold on
                
            end
            title(['Ch',num2str(ch_layout(ch)),';R2:',num2str(R2_Ch_mean,'%1.2f')]); 
            
        end
        
        
    end
    figure(Fe),
    Bar=colorbar('Position', [0.92  0.11  0.01  .84]);
    set(Bar, 'ylim', [-1 1])
    Bar.Ticks = [-1 0 1];
    Bar.TickLabels = {num2str(aclim),num2str(0),num2str(bclim)};
    HighQualityFigs([FeatureName{Fe},'_WeightsVxFixed_6Days(May14-15-21-24-29-31)'])
    
    figure(Fe+NFeatures),
    Bar=colorbar('Position', [0.92  0.11  0.01  .84]);
    set(Bar, 'ylim', [-1 1])
    Bar.Ticks = [-1 0 1];
    Bar.TickLabels = {num2str(aclim),num2str(0),num2str(bclim)};
    HighQualityFigs([FeatureName{Fe},'_WeightsVyFixed_6Days(May14-15-21-24-29-31)'])   
    
    
end


%%



