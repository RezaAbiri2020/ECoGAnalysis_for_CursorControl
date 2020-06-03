

%% This is for calculating Averaged Features per target for imagined blocks
clear all
close all
NCh=128;

%% loading the imagined file or adapt file or fixed file

%load ('..\Data\20190426\GangulyServer\Center-Out\20190426\111618\BCI_Fixed\Data0001.mat')

datadir = uigetdir();
datafiles = dir(fullfile(datadir,'Data*.mat'));

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

for trial=1:length(datafiles)
    %trial=1:NBlocks*8
    %1:length(datafiles)
    %(length(datafiles)-(NBlocks*8-1)):length(datafiles) 
    
    % load data and classifying based on TargetID
    load(fullfile(datadir,datafiles(trial).name))
    
    % Finding the targetID
    TargetID=TrialData.TargetID;
    % Finding the duration of neurofeedback
    idx=find(TrialData.Time>=TrialData.Events(2).Time);
    
    % For different blocks; the extracted neural feature data will be different
    if contains(datadir,'Imagined')
        
        Nbin=length(idx);
        
    elseif contains(datadir,'BCI_Fixed') || contains(datadir,'BCI_CLDA') % one feature value per bin
           
        if TrialData.NeuralSamps(10)<150 % binsize is 100ms
            Nbin=20-5;
            
        elseif TrialData.NeuralSamps(10)>150 % binsize is 200ms
            Nbin=10;
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
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter1).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter1).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==2
        trial_counter2=trial_counter2+1;
        CursorStates(TargetID).Targets(trial_counter2).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter2).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter2).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter2).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter2).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==3
        trial_counter3=trial_counter3+1;
        CursorStates(TargetID).Targets(trial_counter3).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter3).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter3).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter3).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter3).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==4
        trial_counter4=trial_counter4+1;
        CursorStates(TargetID).Targets(trial_counter4).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter4).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter4).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter4).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter4).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==5
        trial_counter5=trial_counter5+1;
        CursorStates(TargetID).Targets(trial_counter5).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter5).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter5).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter5).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter5).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==6
        trial_counter6=trial_counter6+1;
        CursorStates(TargetID).Targets(trial_counter6).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter6).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter6).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter6).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter6).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==7
        trial_counter7=trial_counter7+1;
        CursorStates(TargetID).Targets(trial_counter7).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter7).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter7).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter7).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter7).Trials=TrialData.NeuralSamps;
        
    elseif TargetID==8
        trial_counter8=trial_counter8+1;
        CursorStates(TargetID).Targets(trial_counter8).Trials=TrialData.CursorState(:,idx);
        ECoGs=TrialData.BroadbandData(:,idx);
        ECoG(TargetID).Targets(trial_counter8).Trials=cell2mat(ECoGs');
        SampleTimes(TargetID).Targets(trial_counter8).Trials=TrialData.Time(:,idx)-TrialData.Events(2).Time;
        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
        NeuralFeatures(TargetID).Targets(trial_counter8).Trials=cell2mat(NeuralFeature);
        SampleNumPoints(TargetID).Targets(trial_counter8).Trials=TrialData.NeuralSamps;
        
    end
    
    
end

%% Calculating averaged feature activities (7 features) per target and plotting on brain
%Feature 1: Phase of delta
%Feature 2: delta power
%Feature 3: theta power
%Feature 4: alpha power
%Feature 5: beta  power
%Feature 6: gamma1 power
%Feature 7: gamma2 power

NFeatures=7;
FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};

   
    
%% Plotting averaged features; all features per target in one figure 
for TargetID=1:8 %number of targets
    
    AllFeatures=cat(2,NeuralFeatures(TargetID).Targets(1:length(NeuralFeatures(TargetID).Targets)).Trials);
    MeanFeatures(:,TargetID)=mean(AllFeatures,2);
    ReshapeFeatures=reshape(MeanFeatures(:,TargetID),NCh,NFeatures);
    
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
    
    % I have Coeffs_Brain for this structure for each feature:
    % 1 2 3 ..
    % 17 18....
    % ........128
    
    % I Should do numbering based on following structure for feeding into brain plots
    % codes
    % 113.....128
    %     .......
    % 1 2 3  ...16
    Features_Brain_Mo=zeros(NCh,NFeatures);
    
    for i=1:NFeatures
        A=reshape(Features_Brain(:,i),16,8);
        A=A';
        A=flip(A);
        A=A';
        A=A(:);
        Features_Brain_Mo(:,i)=A;
    end
    
    % plot the features on the brain
    
    load('BRAVO1_lh_pial')
    load('elecs_all')
    
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['averaged features activites (during imagined blocks) on brain: Target ',num2str(TargetID)])
    for i=1:NFeatures
        subplot(2,4,i)
        aclim=-1;
        bclim=1;
        clims = [aclim bclim];
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Features_Brain_Mo(:,i),'lh');
	    %caxis(clims)
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (FeatureName(i))
        %colorbar
    end
    
    %HighQualityFigs(['AveFeaturesTarget',num2str(TargetID),'_Poor'])
    
end

%% Plotting averaged features for all features for all targets in one figure
figure,
set(gcf, 'Position', [100, 100, 2400, 1200]);
suptitle('averaged features activites (during imagined blocks) on brain for all Targets; Sucess rate')

for TargetID=1:8 %number of targets
    
    AllFeatures=cat(2,NeuralFeatures(TargetID).Targets(1:length(NeuralFeatures(TargetID).Targets)).Trials);
    MeanFeatures(:,TargetID)=mean(AllFeatures,2);
    ReshapeFeatures=reshape(MeanFeatures(:,TargetID),NCh,NFeatures);
    
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
    
    % I have Coeffs_Brain for this structure for each feature:
    % 1 2 3 ..
    % 17 18....
    % ........128
    
    % I Should do numbering based on following structure for feeding into brain plots
    % codes
    % 113.....128
    %     .......
    % 1 2 3  ...16
    Features_Brain_Mo=zeros(NCh,NFeatures);
    
    for i=1:NFeatures
        A=reshape(Features_Brain(:,i),16,8);
        A=A';
        A=flip(A);
        A=A';
        A=A(:);
        Features_Brain_Mo(:,i)=A;
    end
    
    % plot the features on the brain
    
    load('BRAVO1_lh_pial')
    load('elecs_all')
    
    
    for i=1:NFeatures
        subplot(8,NFeatures,NFeatures*(TargetID-1)+i)
        aclim=-1;
        bclim=1;
        clims = [aclim bclim];
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Features_Brain_Mo(:,i),'lh');
	    %caxis(clims)
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (['Target: ',num2str(TargetID),'; ',FeatureName(i)])
        %colorbar
    end
    
    %HighQualityFigs(['AveFeaturesTarget',num2str(TargetID),'_Poor'])
    
end



%% Plotting averaged features per feature for all targets in one figure 

for TargetID=1:8 %number of targets
    
    AllFeatures=cat(2,NeuralFeatures(TargetID).Targets(1:length(NeuralFeatures(TargetID).Targets)).Trials);
    MeanFeatures(:,TargetID)=mean(AllFeatures,2);
    ReshapeFeatures=reshape(MeanFeatures(:,TargetID),NCh,NFeatures);
    
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
    
    % I have Coeffs_Brain for this structure for each feature:
    % 1 2 3 ..
    % 17 18....
    % ........128
    
    % I Should do numbering based on following structure for feeding into brain plots
    % codes
    % 113.....128
    %     .......
    % 1 2 3  ...16
    Features_Brain_Mo=zeros(NCh,NFeatures);
    
    for j=1:NFeatures
        A=reshape(Features_Brain(:,j),16,8);
        A=A';
        A=flip(A);
        A=A';
        A=A(:);
        Features_Brain_Mo(:,j)=A;
    end
    
    AllTargets(TargetID).Features=Features_Brain_Mo;
    
end

for i=1:NFeatures
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['averaged features activites (during imagined blocks; Sucess rate) for all targets: ',FeatureName(i)])
    load('BRAVO1_lh_pial')
    load('elecs_all')
    
    for TargetID=1:8
        Features_PerTarget=AllTargets(TargetID).Features;
        % plot the features on the brain
        subplot(2,4,TargetID)
        aclim=-1;
        bclim=1;
        clims = [aclim bclim];
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Features_PerTarget(:,i),'lh');
        %caxis(clims)
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (['Target: ',num2str(TargetID)])
        %colorbar
    end
    
    %HighQualityFigs(['AveFeaturesTarget',num2str(TargetID),'_Poor'])
    
end

%% Plotting averaged features per feature for all targets in one figure 

for TargetID=1:8 %number of targets
    
    AllFeatures=cat(2,NeuralFeatures(TargetID).Targets(1:length(NeuralFeatures(TargetID).Targets)).Trials);
    MeanFeatures(:,TargetID)=mean(AllFeatures,2);
    ReshapeFeatures=reshape(MeanFeatures(:,TargetID),NCh,NFeatures);
    
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
    
    % I have Coeffs_Brain for this structure for each feature:
    % 1 2 3 ..
    % 17 18....
    % ........128
    
    % I Should do numbering based on following structure for feeding into brain plots
    % codes
    % 113.....128
    %     .......
    % 1 2 3  ...16
    Features_Brain_Mo=zeros(NCh,NFeatures);
    
    for j=1:NFeatures
        A=reshape(Features_Brain(:,j),16,8);
        A=A';
        A=flip(A);
        A=A';
        A=A(:);
        Features_Brain_Mo(:,j)=A;
    end
    
    AllTargets(TargetID).Features=Features_Brain_Mo;
    
end

for i=1:NFeatures
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['averaged features activites (during imagined blocks; Sucess rate) for all targets: ',FeatureName(i)])
    load('BRAVO1_lh_pial')
    load('elecs_all')
    
    for TargetID=1:8
        Features_PerTarget=AllTargets(TargetID).Features;
        % plot the features on the brain
        subplot(2,4,TargetID)
        aclim=-1;
        bclim=1;
        clims = [aclim bclim];
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Features_PerTarget(:,i),'lh');
        %caxis(clims)
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (['Target: ',num2str(TargetID)])
        %colorbar
    end
    
    %HighQualityFigs(['AveFeaturesTarget',num2str(TargetID),'_Poor'])
    
end

%% (1) Tuning the values for all featues for all targets on grid to do selection for chs and features...

% calculating all features, their means, STD, for all targets and save in a matrix
for TargetID=1:8 %number of targets
    
    AllFeatures=cat(2,NeuralFeatures(TargetID).Targets(1:length(NeuralFeatures(TargetID).Targets)).Trials);
    MeanFeatures(:,TargetID)=mean(AllFeatures,2);
    ReshapeFeatures=reshape(MeanFeatures(:,TargetID),NCh,NFeatures);
    
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
    
    % I have Coeffs_Brain for this structure for each feature:
    % 1 2 3 ..
    % 17 18....
    % ........128
    
    All_Targets(TargetID).Features=Features_Brain;
    
    
end

% Re-arranging the matrix
for i=1:NFeatures
    
    for TargetID=1:8
        Features_PerTarget=All_Targets(TargetID).Features;
        All_Features(i).Targets(:,TargetID)=Features_PerTarget(:,i);
        
    end
end

% normalizing per channel for each feature
for i=1:NFeatures
    
    for ch=1:NCh
        All_Features_NormPerCh(i).Targets(ch,:)=normalize(All_Features(i).Targets(ch,:),'range',[-1,1]);
    end
end

% normalizing per grid for each feature 
for i=1:NFeatures
    M=All_Features(i).Targets;
    M=M(:);
    M=normalize(M,'range',[-1,1]);
    M=reshape(M,NCh,8);
    All_Features_NormPerGrid(i).Targets=M;
    
end

% ploting normalized per channel for each feature
for i=1:NFeatures
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['contribution of averaged feature activity (Normalized per Ch) from imagined blocks per Ch for all targets: ',FeatureName(i)])
    
    for ch=1:NCh
        
        for TargetID=1:8
            Feature_value=All_Features_NormPerCh(i).Targets(ch,TargetID);
            x=abs(Feature_value)*cos(-(TargetID-1)*2*pi/8);
            y=abs(Feature_value)*sin(-(TargetID-1)*2*pi/8);
            % plot the features on the brain
            subplot(8,16,ch)
            if Feature_value >=0
                plot([0,x],[0,y],'-b')
                %p=compass(x,y);
                % a=annotation('arrow',[0 x],[0 y]);
            else % Feature_value <0
                plot([0,x],[0,y],'-r')
                %p=compass(x,y,'r');
                %a=annotation('arrow',[0 x],[0 y]);
                %a.Color='red';
            end
            ylim([-1 1])
            xlim([-1 1])
            axis off
            %gid off
            hold on
            % title (['Target: ',num2str(TargetID)
            if TargetID==8
                % find all of the lines in the polar plot
                %h = findall(gcf,'type','line');
                % remove the handle for the polar plot line from the array
                %h(h == p) = [];
                % delete all other lines
                %delete(h);
                % find all of the text objects in the polar plot
                %t = findall(gcf,'type','text');
                % delete the text objects
                %delete(t);
                
            end
            
        end
        hold on;
        ezplot('x^2+y^2=1');
        title(['Ch',num2str(ch_layout(ch))]);
        
    end
    
    HighQualityFigs(['AveFeTuningNPerCh',num2str(i),'_Good'])
    
end

%% (2) Tuning the values for all featues for all targets on grid to do selection for chs and features...

% calculating all features, their means, STD, for all targets and save in a matrix
for TargetID=1:8 %number of targets
    
    for trial=1:length(NeuralFeatures(TargetID).Targets)
        Features_PerTgTr=mean(NeuralFeatures(TargetID).Targets(trial).Trials,2);
        
        
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
        
        % I have Coeffs_Brain for this structure for each feature:
        % 1 2 3 ..
        % 17 18....
        % ........128
        
        All_Features(TargetID).Targets(trial).Trials=Features_Brain;
        %All_Features(TargetID).Targets(trial).Trials=abs(Features_Brain);
        
    end
    
end

% calculating the regress/tuning for all targets
% f=B*X
angles=[0:-45:-(360-45)];
X_1=[ones(1,8);sind(angles);cosd(angles)];
X=repmat(X_1,1,length(All_Features(1).Targets));

% each feature on one figure
for i=1:NFeatures
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{i},'; FixedBlocks: Tuning of each Ch [R2>=0.08](the green circle has a unit radius)'])
    
    for ch=1:NCh
        
        f=[];
        for TargetID=1:8
            
            for trial=1:length(All_Features(TargetID).Targets)
                f=[f; All_Features(TargetID).Targets(trial).Trials(ch,i)];
            end
        
        end
        
        [b,bint,r,rint,stats] = regress(f,X');
        if 0
            figure;
            set(gcf, 'Position', [100, 100, 2400, 600]);
            plot(f)
            hold on; plot(b'*X)
            vline(0:8:80)
            xticks([1:1:80])
            xticklabels({'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8',...
                'T1','T2','T3','T4','T5','T6','T7','T8'})
            title(['Tuning of Ch ',num2str(ch_layout(ch)),' for ',FeatureName{i},' with R2:',num2str(stats(1),'%1.2f')]);
            xlabel('Targets')
            ylabel('Feature value')
            HighQualityFigs(['phase of delta_Apr26_Imagined_TuningCh',num2str(ch_layout(ch))])
        
        end
        
        PD=atan2(b(2), b(3)); % radian; to change to degree: /pi*180
        DepthModul=sqrt(b(2)^2+b(3)^2)/abs(b(1));
        R2=stats(1);
        
        if R2>=0.08
            subplot(8,16,ch)
            plot([0 DepthModul*cos(PD)],[0 DepthModul*sin(PD)],'-')
            %compass(DepthModul*cos(PD),DepthModul*sin(PD))
            limit=max(abs(DepthModul*cos(PD)),abs(DepthModul*sin(PD)));
            if limit<1
                limit=1;
            end
            ylim([-limit +limit])
            xlim([-limit +limit])
            %polarplot([0 PD],[0 DepthModul],'-')
            hold on
            ezplot('x^2+y^2=1');
            axis off
            title(['Ch',num2str(ch_layout(ch)),';R2:',num2str(R2,'%1.2f')]);
        else
            subplot(8,16,ch)
            ezplot('x^2+y^2=1');
            axis off
            title(['Ch',num2str(ch_layout(ch))]);
        end
        
        
    end
    
    %HighQualityFigs([FeatureName{i},'_Apr26_Imagined_TuningChs']) 
    %HighQualityFigs([FeatureName{i},'_Apr26_Adapt_TuningChs'])
    HighQualityFigs([FeatureName{i},'_Apr26_Fixed_TuningChs'])
    
end

% all features on one figure
figure,
set(gcf, 'Position', [100, 100, 2400, 1200]);
suptitle('Tuning of each Ch for all features: F1:B F2:R F3:Y F4:G F5:K F6:C F7:M')

for ch=1:NCh
    
    for i=1:NFeatures
        
        f=[];
        for TargetID=1:8
            
            for trial=1:length(All_Features(TargetID).Targets)
                f=[f; All_Features(TargetID).Targets(trial).Trials(ch,i)];
            end
            
        end
        
        [b,bint,r,rint,stats] = regress(f,X');
        PD(i)=atan2(b(2), b(3)); % radian; to change to degree: /pi*180
        DepthModul(i)=sqrt(b(2)^2+b(3)^2)/abs(b(1));
        R2(i)=stats(1);
        
    end
    
    subplot(8,16,ch)
    hold on
    u=[zeros(7,1), (DepthModul.*cos(PD))'];
    v=[zeros(7,1), (DepthModul.*sin(PD))'];
    plot(u(1,:),v(1,:),'-b',u(2,:),v(2,:),'-r',u(3,:),v(3,:),'-y',u(4,:),v(4,:),'-g',...
        u(5,:),v(5,:),'-k',u(6,:),v(6,:),'c',u(7,:),v(7,:),'m-')
    
    limit=max(max(abs(u(:,2)),abs(v(:,2))));
    if limit<1
        limit=1;
    end
    ylim([-limit +limit])
    xlim([-limit +limit])
    hold on
    ezplot('x^2+y^2=1');
    axis off
    title(['Ch',num2str(ch_layout(ch)),';R2:',num2str(mean(R2),'%1.2f')]);
    
    
end

 %HighQualityFigs('TuningChsForAllFeatures_Poor')
 HighQualityFigs('TuningImagined(Apr26-5hz5hz_Good)')
 HighQualityFigs('TuningAdapt_AllBlocks(Apr26-5hz5hz_Good)')
 HighQualityFigs('TuningFixed(Apr26-5hz5hz_Good)')  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx1=TrialData.Time-TrialData.TrialStartTime
idx2=TrialData.Time-TrialData.Events(1).Time % instructed delay
idx3=TrialData.Time-TrialData.Events(2).Time % Reach Target
idx4=TrialData.TrialEndTime-TrialData.Time