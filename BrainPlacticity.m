

%% Tuning analysis for the first 2 seconds of data in fixed blocks

% Performing compass plots using R2 and PD
% Performing compass plots using DM and PD

% the program can be run for imagined or adapt blocks as well
clear all
close all


%% including the selected Fixed sessions for each day for analysis through

Days={'20190401','20190403','20190412',...
    '20190423','20190426','20190429',...
    '20190501','20190507','20190510','20190514','20190515',...
    '20190521','20190524','20190529','20190531','20190604',...
    '20190607','20190610'};
% corresponding sessions for mentioned days

Sessions(1).Days={'111316', '111942', '112738', '133838', '135810', '141638'};
Sessions(2).Days={'105510', '111944', '132634', '135357'};
Sessions(3).Days={'111554'};
Sessions(4).Days={'142332', '144207', '144915','150325', '151756'};
Sessions(5).Days={'111618', '115111', '141302', '144334'};
Sessions(6).Days={'130907', '134613', '135428', '140100'};
Sessions(7).Days={'105345', '112049','135512'};
Sessions(8).Days={'111526', '113944', '140915'};
Sessions(9).Days={'111234', '140202', '141900'};
Sessions(10).Days={'135927', '140808', '141731', '142704'};
Sessions(11).Days={'112302', '113447', '135831', '141134', '141859'};
Sessions(12).Days={'135943', '141438'};
Sessions(13).Days={'111731', '112353','134653', '135957'};
Sessions(14).Days={'111528', '114050'};
Sessions(15).Days={'105048', '110703', '112444', '133517', '140204', '141319'};
Sessions(16).Days={'112454', '114636', '141420', '143109'};
Sessions(17).Days={'105615', '135050', '140828'};
Sessions(18).Days={'105809', '110746', '132416', '133734', '135334'};


for day=1:length(Days)
    direc=['../Data\',char(Days(day)),'\GangulyServer\Center-Out\',char(Days(day)),'\'];
    Ses=Sessions(day).Days;
    Bravo1(day).Days=cell(1,length(Ses));
    
    for S=1:length(Ses)
        Bravo1(day).Days(1,S)={[direc,char(Ses(1,S)),'\BCI_Fixed']};
    end
    
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
        
        if strcmp(class(TrialData.KalmanFilter),'cell')
            
            Coeff_Features=TrialData.KalmanFilter{1,end}.C;
            
        elseif strcmp(class(TrialData.KalmanFilter),'struct')
            
            Coeff_Features=TrialData.KalmanFilter.C;
            
        end
  
        % third column for vx and forth column for vy
        %plot(Coeff_Features(:,3))
        %plot(Coeff_Features(:,4))
        
        NCh=TrialData.Params.NumChannels;
        
        % for vx
        Vx_Coeffs=reshape(Coeff_Features(:,3),NCh,NFeatures);
        % for vy
        Vy_Coeffs=reshape(Coeff_Features(:,4),NCh,NFeatures);
        
        % organizing the weight for plotting
        ch_layout = [
            96    84    76    95    70    82    77    87    74    93    66    89    86    94    91    79
            92    65    85    83    68    75    78    81    72    69    88    71    80    73    90    67
            62    37    56    48    43    44    60    33    49    64    58    59    63    61    51    34
            45    53    55    52    35    57    38    50    54    39    47    42    36    40    46    41
            19    2    10    21    30    23    17    28    18    1    8    15    32    27    9    3
            24    13    6    4    7    16    22    5    20    14    11    12    29    26    31    25
            124    126    128    119    110    113    111    122    117    125    112    98    104    116    103    106
            102    109    99    101    121    127    105    120    107    123    118    114    108    115    100    97];
        
        
        ch_layout_1=ch_layout';
        ch_layout_2=ch_layout_1(:);
        Vx_Coeffs_Brain=zeros(NCh,NFeatures);
        Vy_Coeffs_Brain=zeros(NCh,NFeatures);
        
        for j=1:NCh
            Vx_Coeffs_Brain(j,:)=Vx_Coeffs(ch_layout_2(j),:);
            Vy_Coeffs_Brain(j,:)=Vy_Coeffs(ch_layout_2(j),:);
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
        
        Hiting_Nu=0;
        Trials_Nu=0;
        
        for trial=[1:length(datafiles)]
            %trial=1:NBlocks*8
            %1:length(datafiles)
            %(length(datafiles)-(NBlocks*8-1)):length(datafiles)
            % load data and classifying based on TargetID
            load(fullfile(datadir,datafiles(trial).name))
            
            if TrialData.ErrorID==0 % he hit the target
                Hiting_Nu=Hiting_Nu+1;
            end 
            
            if ~isempty(TrialData.Params.BadChannels)
                
                fprintf(['Bad Chs for \n day:',char(Days(day)),'\n Session:',char(Sessions(day).Days(1,session)),'\n']);
                
            end
            
            if length(TrialData.Events)==2 % it means the cursor start from center to hit the targets
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
                Trials_Nu=Trials_Nu+1;
                % Finding the targetID
                TargetID=TrialData.TargetID;
                % Finding the duration of neurofeedback
                             
                Nbin_Hit=length(idx);
                TimeForAnalysis=2; % second
                
                if TrialData.Params.UpdateRate==10 %TrialData.NeuralSamps(5)<150 % binsize is 100ms
                    Nbin=(10*TimeForAnalysis);
                    
                elseif TrialData.Params.UpdateRate==5 %TrialData.NeuralSamps(5)>150 % binsize is 200ms
                    Nbin=(5*TimeForAnalysis);
                end
                
                
                % Extracting data: each row will be for a target with corresponding
                % trials inside that
                if TargetID==1
                    trial_counter1=trial_counter1+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter1).Trials=cell2mat(NeuralFeature);
                    
                elseif TargetID==2
                    trial_counter2=trial_counter2+1;
                                       
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter2).Trials=cell2mat(NeuralFeature);
                    
                    
                elseif TargetID==3
                    trial_counter3=trial_counter3+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter3).Trials=cell2mat(NeuralFeature);
                    
                elseif TargetID==4
                    trial_counter4=trial_counter4+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter4).Trials=cell2mat(NeuralFeature);
                    
                elseif TargetID==5
                    trial_counter5=trial_counter5+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter5).Trials=cell2mat(NeuralFeature);
                    
                elseif TargetID==6
                    trial_counter6=trial_counter6+1;
                                       
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter6).Trials=cell2mat(NeuralFeature);
                   
                elseif TargetID==7
                    trial_counter7=trial_counter7+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter7).Trials=cell2mat(NeuralFeature);
                    
                elseif TargetID==8
                    trial_counter8=trial_counter8+1;
                    
                    if Nbin_Hit>=Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin));
                    elseif Nbin_Hit<Nbin
                        NeuralFeature=TrialData.NeuralFeatures(:,idx(1):idx(Nbin_Hit));
                        NeuralFeature(1,Nbin_Hit+1:Nbin)={NaN(size(NeuralFeature{Nbin_Hit}))};
                    end
                    
                    NeuralFeatures(TargetID).Targets(trial_counter8).Trials=cell2mat(NeuralFeature);
                    
                end
                
            end
        end
        
        % Fining the values for Tuning ..
        
        % calculating all features, their means, STD, for all targets and save in a matrix
        for TargetID=1:length(NeuralFeatures) %number of targets
            
            for trial=1:length(NeuralFeatures(TargetID).Targets)
                Features_PerTgTr=nanmean(NeuralFeatures(TargetID).Targets(trial).Trials,2);
                
                ReshapeFeatures=reshape(Features_PerTgTr,NCh,NFeatures);
                
                Features_Brain=zeros(NCh,NFeatures);
                
                for j=1:NCh
                    Features_Brain(j,:)=ReshapeFeatures(ch_layout_2(j),:);
                end
                
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
                DepthModul(ch)=sqrt(b(2)^2+b(3)^2);%/abs(b(1));
                R2(ch)=stats(1);
                
            end
            AllTunings(day).Days(session).Sessions(Fe).Features.PD=PD;
            AllTunings(day).Days(session).Sessions(Fe).Features.DepthModul=DepthModul;
            AllTunings(day).Days(session).Sessions(Fe).Features.R2=R2;
            
            
        end
        
        clear Vx_Coeffs_Brain Vy_Coeffs_Brain;
        clear NeuralFeatures All_Features
        
    SRate=Hiting_Nu/Trials_Nu;
    AllTunings(day).Days(session).SessionsSuccessRate=SRate;
    
    end

    
end

%% Plotting compass plots for magnitude R2 and its angle PD

% plotting across days for individual Ch by considering the features. If
% there are 7 features, then plot. Otherwise ignore that session or day.
%Feature 1: Phase of delta
%Feature 2: delta power
%Feature 3: theta power
%Feature 4: alpha power
%Feature 5: beta  power
%Feature 6: gamma1 power
%Feature 7: gamma2 power

FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};

for Fe=1:7
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{Fe},'; compass plots for magnitude R2 and its angle PD(tuning) inside each Ch across days'])
    
    for ch=1:128
        
        for day=1:length(AllTunings)
            
            for session=1:length(AllTunings(day).Days)
                
                if length(AllTunings(day).Days(session).Sessions)==7
                    
                    Magnitude=AllTunings(day).Days(session).Sessions(Fe).Features.R2(ch);
                    CalAngle=AllTunings(day).Days(session).Sessions(Fe).Features.PD(ch);
                 
                    subplot(8,16,ch)
                    Clr=(day*50+session*8)/1000;
                    patchline([0 Magnitude*cos(CalAngle)],[0 Magnitude*sin(CalAngle)],...
                        'linewidth',.5,'edgecolor','b','facecolor','b','FaceAlpha',Clr,'EdgeAlpha',Clr);
                end
            end
        end
        
        Lim=gca;
        maxlim = max(abs([Lim.XLim Lim.YLim]));
        xlim([-maxlim maxlim]);
        ylim([-maxlim maxlim]);
        %xlim([-.5 .5]);
        %ylim([-.5 .5]);
        %set(gca,'xtick',[])
        %set(gca,'ytick',[])
        box on
        title(['Ch',num2str(ch_layout_2(ch))]);
        
        
    end
    HighQualityFigs([FeatureName{Fe},'_ChangeOfR2DPTuning_AcrossDays'])
end

%% Plotting compass plots for magnitude Depth of Modulation and its angle PD

% plotting across days for individual Ch by considering the features. If
% there are 7 features, then plot. Otherwise ignore that session or day.
%Feature 1: Phase of delta
%Feature 2: delta power
%Feature 3: theta power
%Feature 4: alpha power
%Feature 5: beta  power
%Feature 6: gamma1 power
%Feature 7: gamma2 power

FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};

for Fe=1:7
    figure,
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{Fe},'; compass plots for magnitude Depth of Modulation and its angle PD(tuning) inside each Ch across days'])
    
    for ch=1:128
        
        for day=1:length(AllTunings)
            
            for session=1:length(AllTunings(day).Days)
                
                if length(AllTunings(day).Days(session).Sessions)==7
                    
                    Magnitude=AllTunings(day).Days(session).Sessions(Fe).Features.DepthModul(ch);
                    CalAngle=AllTunings(day).Days(session).Sessions(Fe).Features.PD(ch);
                 
                    subplot(8,16,ch)
                    Clr=(day*50+session*8)/1000;
                    patchline([0 Magnitude*cos(CalAngle)],[0 Magnitude*sin(CalAngle)],...
                        'linewidth',.5,'edgecolor','b','facecolor','b','FaceAlpha',Clr,'EdgeAlpha',Clr);
                end
            end
        end
        
        Lim=gca;
        maxlim = max(abs([Lim.XLim Lim.YLim]));
        xlim([-maxlim maxlim]);
        ylim([-maxlim maxlim]);
        %set(gca,'xtick',[])
        %set(gca,'ytick',[])
        box on
        title(['Ch',num2str(ch_layout_2(ch))]);
        
        
    end
    HighQualityFigs([FeatureName{Fe},'_ChangeOfDMDPTuning_AcrossDays'])
end


