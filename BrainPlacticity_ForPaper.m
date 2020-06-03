

%% Tuning analysis for the first 2 seconds of data in fixed blocks
% Plotting R2, Depth of Modulation and Preferred Direction distribution for all days and for specific chs
%  weights are also recored for all days and for specific chs

% the program cab be run for imagined or adapt blocks as well
clear all
close all



%% including the selected Fixed sessions for each day for analysis through

Days={'20190501', '20190507','20190510','20190514','20190515','20190521'...
    ,'20190524','20190529','20190531','20190604',...
    '20190607','20190610','20190618','20190621','20190626'};
% corresponding sessions for mentioned days
Sessions(1).Days={'105345','112049','135512'};
Sessions(2).Days={'111526','113944'};
Sessions(3).Days={'141900'};
Sessions(4).Days={'135927', '140808', '141731', '142704'};
Sessions(5).Days={'112302', '113447'};
Sessions(6).Days= {'135943','141438'};
Sessions(7).Days={'111731', '112353','134653', '135957'};
Sessions(8).Days={'111528', '114050'};
Sessions(9).Days={'105048', '110703', '112444', '133517', '140204', '141319'};
Sessions(10).Days={'112454', '114636', '143109'};
Sessions(11).Days={'105615', '135050', '140828'};
Sessions(12).Days={'105809', '110746'};
Sessions(13).Days= {'133933','135944','143103'};
Sessions(14).Days={'110942','113715','134143','141129'};
Sessions(15).Days={'105539','110229','110553','111358','112324'};

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
                PD=atan2(b(2), b(3)); % radian; to change to degree: /pi*180
                DepthModul=sqrt(b(2)^2+b(3)^2)/abs(b(1));
                R2=stats(1);
                
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.f=f;
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.X=X;
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.b=b;
                                
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.PD=PD;
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.DepthModul=DepthModul;
                AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.R2=R2;
                
            end
   
        end
        
        %AllTunings(day).Days(session).Sessions.Counter=All_Counters;
        
        clear Vx_Coeffs_Brain Vy_Coeffs_Brain;
        clear NeuralFeatures All_Features
        
    SRate=Hiting_Nu/Trials_Nu;
    AllTunings(day).Days(session).SessionsSuccessRate=SRate;
    AllTunings(day).Days(session).CounterTrialsPerTarget=All_Counters;
                
    
    end

    
end

%% ploting tuning for a specific Fe & Ch

% selected Chs on grid based on (BrainPlacticity_V1_ForPaper_V0) for delta band 
%Fe=2;
%ch=62; % ch 40 inside blackrock
%ch=95; % ch 103 inside blackrock 
%ch=61; % ch 36 inside blackrock


% selected Chs on grid based on (BrainPlacticity_V1_ForPaper_V0) for High gamma band 
Fe=7;
ch=96; % ch 106 inside blackrock
%ch=128; % ch 97 inside blackrock 
%ch=80; % ch 3 inside blackrock
%ch=94; % ch 26 inside blackrock

% for saving the fig
ChName=['Ch36','(',num2str(ch),')'];

clear F_train
for day=1:length(Bravo1)
    
    for session=1:length(Bravo1(day).Days)
        
        X=AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.X;
        f=AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.f;
        b=AllTunings(day).Days(session).Sessions(Fe).Features(ch).Channels.b;
        Counter=AllTunings(day).Days(session).CounterTrialsPerTarget;
        
        if Counter(5)==0 && Counter(6)==0 && Counter(7)==0 && Counter(8)==0
            
            F1=NaN(max(Counter),4);
            if Counter(1)~=0
                F1(1:Counter(1),1)=f(1:Counter(1));
            end
            
            for kk=2:4
                if Counter(kk)~=0
                    F1(1:Counter(kk),kk)=f(sum(Counter(1:(kk-1)))+1:sum(Counter(1:kk)));
                end
            end
            
            angles_day=[0;-90;-180;-270];
            All_angles=[0:-45:-315]';
            
            vq=[];
            for k=1:max(Counter)
            vq(k,:)= interp1(angles_day,F1(k,:),All_angles,'pchip');
            end
            
            F_train(day).Days(session).Sessions=vq;
        
        else
            
            F1=NaN(max(Counter),8);
            if Counter(1)~=0
            F1(1:Counter(1),1)=f(1:Counter(1));
            end 
            
            for kk=2:8
                if Counter(kk)~=0
                    F1(1:Counter(kk),kk)=f(sum(Counter(1:(kk-1)))+1:sum(Counter(1:kk)));
                end
            end
            
            All_angles=[0:-45:-315]';
            
            vq=[];
            for k=1:max(Counter)
                vq(k,:)= interp1(All_angles,F1(k,:),All_angles,'pchip');
            end
            
            F_train(day).Days(session).Sessions=vq; 
             
        end
        
        
    end 
end 

%% plot tuning for all days with error bar for sessions

clear F_predictDay_NonSmooth F_predictDay_Smooth

R2_Ch=[];PD_ch=[];
for day=1:length(Bravo1)
    
    F_PerDay=[];
    X=[]; X_day=[];
    
    for session=1:length(Bravo1(day).Days)
        F_PerDay=[F_PerDay; F_train(day).Days(session).Sessions];
        
    end
    F_PerDay=F_PerDay';
    F_PerDay=F_PerDay(:);
    X=[ones(1,8);sind(0:-45:-315);cosd(0:-45:-315)];
    X_day=repmat(X,1,length(F_PerDay)/8);
    [b,bint,r,rint,stats] = regress(F_PerDay,X_day');
    R2_Ch(day)=stats(1);
    PD_Ch(day)=atan2(b(2), b(3))/pi*180; % radian; to change to degree: /pi*180
    F_predictDay_NonSmooth(day).Days=b'*X;
    X_Smooth=[ones(1,360);sind(0:-1:-359);cosd(0:-1:-359)];
    F_predictDay_Smooth(day).Days=b'*X_Smooth;
    
end

% plot: smooth cosine tuning for 15days on (one figure) with blue color for manuscript 

figure;
set(gcf, 'Position', [1500, 800, 500, 400]);
ColorCodes=brewermap(15,'Blues');
for day=1:15
    plot(1:7/360:(8-7/360),F_predictDay_Smooth(day).Days,'linewidth',1.5,'color',ColorCodes(day,:))
    hold on
end 
set(gca,'FontSize',16)
xlim([.5,8.5])
%ylim([-1,2])
xticks([1,8])
xticklabels({['0',char(176)],['315',char(176)]})
xlabel('Direction','FontSize', 14')
%ylabel('Delta Power (STD)','FontSize', 14')
ylabel('High Gamma Power (STD)','FontSize', 14')
%HighQualityFigs(['TuningsFixed_Delta_',ChName,'_AllDays'])
%HighQualityFigs(['TuningsFixed_HighGamma_',ChName,'_AllDays'])


%% plot "F_predictDay" per day with error bar for "F_train" per day (for selected ch; 15 figures)

for day=2%1:length(F_train)
    R2_Means=[];
    R2_errors=[];
    F_values=[];
    for session=1:length(F_train(day).Days)
        F_values=[F_values; F_train(day).Days(session).Sessions];
    end
    R2_Means=mean(F_values,1);
    R2_errors=std(F_values,1);
    figure;
    set(gcf, 'Position', [1500, 800, 500, 400]);
    errorbar(R2_Means,R2_errors,'rs','MarkerSize',10,...
        'MarkerEdgeColor','red')%,'MarkerFaceColor','red')
    
    hold on;
    %plot a smooth tuning cosine curve
    plot(1:7/360:(8-7/360),F_predictDay_Smooth(day).Days,'-b','linewidth',1.5)
    set(gca,'FontSize',16)
    xlim([.5,8.5])
    %ylim([-2.5,2.5])
    xticks([1,8])
    xticklabels({['0',char(176)],['315',char(176)]})
    xlabel('Direction','FontSize', 14')
    %ylabel('Delta Power (STD)','FontSize', 14')
    ylabel('High Gamma Power (STD)','FontSize', 14')
    %HighQualityFigs(['TuningsFixed_Delta_',ChName,'_Day',num2str(day)])
    %HighQualityFigs(['TuningsFixed_HighGamma_',ChName,'_Day',num2str(day)])
    
end


%% ploting R2 across days for selected ch
figure;
set(gcf, 'Position', [1500, 800, 500, 400]);
plot(R2_Ch,'or','linewidth',2)
hold on
plot(R2_Ch,'-b','linewidth',2)
set(gca,'FontSize',20)
xlim([0,16])
%ylim([-.5,2.5])
xticks([1,5,10,15])
xticklabels({'1','5','10','15'})
xlabel('Day','FontSize', 20')
ylabel('R^2 (STD)','FontSize', 20')
HighQualityFigs(['R2ChangeAcrossDays_Delta_',ChName])
%HighQualityFigs(['R2ChangeAcrossDays_HighGamma_',ChName])


%% ploting PD across days for selected ch
figure;
set(gcf, 'Position', [1500, 800, 500, 400]);
plot(PD_Ch,'or','linewidth',2)
hold on
plot(PD_Ch,'-b','linewidth',2)
set(gca,'FontSize',20)
xlim([0,16])
%ylim([-.5,2.5])
xticks([1,5,10,15])
xticklabels({'1','5','10','15'})
xlabel('Day','FontSize', 20')
ylabel(['Preferred Direction (',char(176),')'],'FontSize', 20')
HighQualityFigs(['PDChangeAcrossDays_Delta_',ChName])
%HighQualityFigs(['PDChangeAcrossDays_HighGamma_',ChName])

