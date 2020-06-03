

%% Plotting tuning and weights on the grid for only fixed blocks (For both center reset on/off/discrete/continuous control sessions).

clear all
close all
NCh=128;

%% Loading the related data set...

datadir = uigetdir();
datafiles = dir(fullfile(datadir,'Data*.mat'));


%% Extracting the variables values for all targets

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

NBlocks=3;

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
            
            if TrialData.NeuralSamps(10)<150 % binsize is 100ms
                Nbin=(10*TimeForAnalysis);
                
            elseif TrialData.NeuralSamps(10)>150 % binsize is 200ms
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

%% Tuning the values for all features for all targets on grid

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

%% Plotting:

figure,
set(gcf, 'Position', [100, 100, 2400, 1200]);
suptitle('May21-Fixed141438: Tuning of each Ch for all features: F1:Blue F2:Red F3:Yellow F4:Green F5:Black F6:Cyan F7:Magenta')

for ch=1:NCh
    
    for i=1:NFeatures
        
        f=[];
        for TargetID=1:length(NeuralFeatures)
            
            for trial=1:length(All_Features(TargetID).Targets)
                f=[f; All_Features(TargetID).Targets(trial).Trials(ch,i)];
            end
            
        end
        
        [b,bint,r,rint,stats] = regress(f,X');
        PD(i)=atan2(b(2), b(3)); % radian; to change to degree: /pi*180
        %DepthModul=sqrt(b(2)^2+b(3)^2)/abs(b(1));
        DepthModul(i)=1;
        R2(i)=stats(1);
        
    end
    subplot(8,16,ch)
    hold on
    u=[zeros(NFeatures,1), (DepthModul.*cos(PD))'];
    v=[zeros(NFeatures,1), (DepthModul.*sin(PD))'];
    plot(u(1,:),v(1,:),'-b',u(2,:),v(2,:),'-r',u(3,:),v(3,:),'-y',u(4,:),v(4,:),'-g',...
        u(5,:),v(5,:),'-k',u(6,:),v(6,:),'c',u(7,:),v(7,:),'m-','linewidth',1.5)
    
    %     limit=max(max(abs(u(:,2)),abs(v(:,2))));
    %     if limit<1
    %     limit=1;
    %     end
    limit=1;
    ylim([-limit +limit])
    xlim([-limit +limit])
    hold on
    ezplot('x^2+y^2=1');
    axis off
    title(['Ch',num2str(ch_layout(ch)),';R2:',num2str(mean(R2),'%1.2f')]);
    
end


HighQualityFigs('TuningAllChsAllFs_May21_Fixed141438')


