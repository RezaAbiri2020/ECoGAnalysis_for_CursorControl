
% see the following for a breif description of this code; Karunesh asked to have similar gap for days forget about trials x axis
% also he asked to caclulate the PnP results for block by block on each
% sesssion

%% check if kalman weigths are converging across days and months for three features setup

% then;
%Section 1:
% plot KF convergence perofrmance with continuation of PnP; their performance; bitrate

%Section 2:
% plot KF convergence performance with continuation of Reset KF test with neural maps


%% for section 1
clear all
close all
clc


%% experiment info for long-term CLDA for section 1; three features: delta, beta, high gamma
expts = [];
 
% Long-Term CLDA in Center-Out using three features: delta, beta, high gamma 
% Testing the performance in Center-Out or Radial tasks, as well
if 0 % for section 1; these days are excluded following Karunesh comments
    
    
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
    
    expts(end+1).yymmdd = '20191122';
    expts(end).CLDA_hhmmss = {'112414','141701'};
    expts(end).perf_CenterOut_hhmmss = {'112414','141701','142524','143443'};
    expts(end).perf_RadialTask_hhmmss = {'113743','114003','114800'};
    
    % from this day using the new ModularBCI codes desgined by Daniel
    expts(end+1).yymmdd = '20191127';
    expts(end).CLDA_hhmmss = {'142835'};
    expts(end).perf_CenterOut_hhmmss = {'142835'};
    expts(end).perf_RadialTask_hhmmss = {'144302','144516','144837','145055','145243','145903'};
end


% RadialTypingMultiClick for the following day
expts(end+1).yymmdd = '20191203';
expts(end).CLDA_hhmmss = {'135420'}; %exclude one session for consistency
% with others in plotting {'134453','135420'}; based on Karunesh comments,
expts(end).perf_CenterOut_hhmmss = {'135420'};%{'134453','135420'};
expts(end).perf_RadialTask_hhmmss = {'141002'};

% RadialTypingMultiClick for the following day
expts(end+1).yymmdd = '20191206';
expts(end).CLDA_hhmmss = {'110809'};
% CLDA blocks without fixed blocks: '140949','142301','142643'
expts(end).perf_CenterOut_hhmmss = {'110809'};
expts(end).perf_RadialTask_hhmmss = {'112558','113508','114445'};

% RadialTyping for the following day
expts(end+1).yymmdd = '20191217';
expts(end).CLDA_hhmmss = {'134525'};
expts(end).perf_CenterOut_hhmmss = {'134525'};
expts(end).perf_RadialTask_hhmmss = {'145444','145838','150222'};

%% Section 1-1: PnP experiment info after previous long-term CLDA (for sections 1-2, 1-3)

expts(end+1).yymmdd = '20200115';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'110726','111243'};

expts(end+1).yymmdd = '20200117';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'112153','113322','114206'};

expts(end+1).yymmdd = '20200124';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'113251','114404','115030','142238'};

expts(end+1).yymmdd = '20200127';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'142333','142707'};

expts(end+1).yymmdd = '20200131';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'111305'};


%% Section 1-2:calculate and plot only the convergence of KF in Long-Term CLDA for three features set for section 1 for PnP

% going through CLDA files and find the angle compared to final reference
Angle_x_AcrossSessions = [];
Angle_y_AcrossSessions = [];
Angle_c_AcrossSessions = [];

Angle_x_WithinSessions = [];
Angle_y_WithinSessions = [];
Angle_c_WithinSessions = [];
Sessions={};

Angle_x_AcrossTrials = [];
Angle_y_AcrossTrials = [];
Angle_c_AcrossTrials = [];

% finding the last trial as reference
Expt_Ref=expts(3); % the last day of CLDA
datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',Expt_Ref.yymmdd,'GangulyServer',...
    Expt_Ref.yymmdd,'CenterOut',Expt_Ref.CLDA_hhmmss{1,length(Expt_Ref.CLDA_hhmmss)},'BCI_CLDA');
files = dir(fullfile(datadir,'Data*.mat'));
load(fullfile(datadir,files(end).name));
C_Ref=TrialData.KalmanFilter{1,1}.C;

% calculate the angles relative to refs
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    if isempty(expt.CLDA_hhmmss)
        Sessions(end+1)={yymmdd};
        DayTrials(i)=DayTrials(i-1)+16; % add 16 trials for plotting purpose
        Angle_x_AcrossTrials= [Angle_x_AcrossTrials;zeros(16,1)];
        Angle_y_AcrossTrials= [Angle_y_AcrossTrials;zeros(16,1)];
        Angle_c_AcrossTrials= [Angle_c_AcrossTrials;zeros(16,1)];
        
    else
        for j=1:length(expt.CLDA_hhmmss)
            Sessions(end+1)={[yymmdd,'-',expt.CLDA_hhmmss{1,j}]};
            
            fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
            
            % go through datafiles in CLDA blocks
                datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,'GangulyServer',...
                    yymmdd,'CenterOut',expt.CLDA_hhmmss{1,j},'BCI_CLDA');
                           
            
            files = dir(fullfile(datadir,'Data*.mat'));
            load(fullfile(datadir,files(1).name));
            C_Start=TrialData.KalmanFilter{1,1}.C;
            
            load(fullfile(datadir,files(end).name));
            C_End=TrialData.KalmanFilter{1,1}.C;
            
            Angle_x_AcrossSessions= [Angle_x_AcrossSessions;subspace(C_End(:,3),C_Ref(:,3))*180/pi];
            Angle_y_AcrossSessions= [Angle_y_AcrossSessions;subspace(C_End(:,4),C_Ref(:,4))*180/pi];
            Angle_c_AcrossSessions= [Angle_c_AcrossSessions;subspace(C_End(:,5),C_Ref(:,5))*180/pi];
            
            Angle_x_WithinSessions= [Angle_x_WithinSessions;subspace(C_End(:,3),C_Start(:,3))*180/pi];
            Angle_y_WithinSessions= [Angle_y_WithinSessions;subspace(C_End(:,4),C_Start(:,4))*180/pi];
            Angle_c_WithinSessions= [Angle_c_WithinSessions;subspace(C_End(:,5),C_Start(:,5))*180/pi];
            
            for k=1:length(files)
                load(fullfile(datadir,files(k).name));
                C_New=TrialData.KalmanFilter{1,1}.C;
                Angle_x_AcrossTrials= [Angle_x_AcrossTrials;subspace(C_New(:,3),C_Ref(:,3))*180/pi];
                Angle_y_AcrossTrials= [Angle_y_AcrossTrials;subspace(C_New(:,4),C_Ref(:,4))*180/pi];
                Angle_c_AcrossTrials= [Angle_c_AcrossTrials;subspace(C_New(:,5),C_Ref(:,5))*180/pi];
                
            end
            
        end
        
        DayTrials(i)=length(Angle_x_AcrossTrials);
    end
    
    
end

% Across trials
figure('position',[400 400 1000 100]);
box off
plot(Angle_x_AcrossTrials,'linewidth',2)
hold on
plot(Angle_y_AcrossTrials,'linewidth',2)
hold on;
ylim([-1 30])
xlim([0 length(Angle_x_AcrossTrials)])
yticks([0,30])
xticks('')
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
vline(DayTrials)
set(gca,'FontSize',20)
%xlabel('Trials')
%title('KF weights convergence across all trials for three features setup')
box off
%HighQualityFigsSVG('Fig_R2_1')


%% Section 1-3: plot KF convergence performance of Long-term CLDA (performance done in different tasks)with continuation of PnP performance in the Radial Tasks

% eliminated for this version; see previous version 1 for codes

%% Section 1-4: plot KF convergence performance of Long-term CLDA in Center-Out with continuation of PnP performance in RadialTask
% Calculate the fitts rates, performance (success rates, time to target) corresponding to CLDA blocks and also PnP tests

Day_Ref='20190521';

% for long-term CLDA & PnP
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    PerfsDays(i).TestDate=yymmdd;
    
    NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
        [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);
    PerfsDays(i).TestDayCounter=NumDays;
    
    % going through the fixed blocks
    if isempty(expt.CLDA_hhmmss)
        
        PerfsDays(i).PnP(1).SessionName='Plug&Play';
        FixedType='RadialTask';
        PerfsDays(i).FixedType=FixedType;
        FixedBlockCounter=0;
        for jj=1:length(expt.perf_RadialTask_hhmmss)            
                FixedBlockCounter=FixedBlockCounter+1;
                %calculate the type of fixed block, fitts rate & performance & clicker
                %on/off
                [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_RadialTask_hhmmss{1,jj},FixedType);
                % save into structure:
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).Clicker=Clicker;     
        end    
        
    else 
        % for each day we only have one type of fixed block; so:
            FixedType='Center-Out';
            PerfsDays(i).FixedType=FixedType;
        
        % go througth the number of CLDA sessions
        for j=1:length(expt.CLDA_hhmmss)
            fprintf('CLDA:')
            fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
            PerfsDays(i).CLDA(j).SessionNameCLDA=expt.CLDA_hhmmss{1,j};
            
            if j<length(expt.CLDA_hhmmss)
                % go througth the number of related fixed sessions
                FixedBlockCounter=0;
                for jj=1:length(expt.perf_CenterOut_hhmmss)
                    
                    if str2num(expt.perf_CenterOut_hhmmss{1,jj})>= str2num(expt.CLDA_hhmmss{1,j}) && ...
                            str2num(expt.perf_CenterOut_hhmmss{1,jj})<str2num(expt.CLDA_hhmmss{1,j+1})
                       
                        FixedBlockCounter=FixedBlockCounter+1;
                        %calculate the type of fixed block, fitts rate & performance & clicker
                        %on/off
                        [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                        % save into structure:
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
                        
                    end
                    
                end
                
            elseif j==length(expt.CLDA_hhmmss)
                
                FixedBlockCounter=0;
                for jj=1:length(expt.perf_CenterOut_hhmmss)
                    
                    if str2num(expt.perf_CenterOut_hhmmss{1,jj})>= str2num(expt.CLDA_hhmmss{1,j})
                        
                        FixedBlockCounter=FixedBlockCounter+1;
                        %calculate the type of fixed block, fitts rate & performance & clicker
                        %on/off
                        [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                        % save into structure:
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
                        
                    end
                    
                end
                
            end
            
        end
    end
    
end

% plotting the results
Perf_Days={};

for i=1:length(PerfsDays)
    
    Perf_Days(i)={num2str(PerfsDays(i).TestDayCounter)};
    counter_CLDA=0;
    
    if isempty(PerfsDays(i).CLDA)% it is PnP
        
        for k=1:length(PerfsDays(i).PnP(1).FixedBlocks)
            counter_CLDA=counter_CLDA+1;
            Perf_Plot(i).SuccessRate(counter_CLDA)=PerfsDays(i).PnP(1).FixedBlocks(k).SuccessRate;
            Perf_Plot(i).TravelTime(counter_CLDA)=PerfsDays(i).PnP(1).FixedBlocks(k).TravelTime;
            Perf_Plot(i).BitRate(counter_CLDA)=PerfsDays(i).PnP(1).FixedBlocks(k).FittsRate;
            Perf_Plot(i).Clicker(counter_CLDA)={PerfsDays(i).PnP(1).FixedBlocks(k).Clicker};
            
        end
        
    else % it is CLDA
        
        for j=1:length(PerfsDays(i).CLDA)
            
            for k=1:length(PerfsDays(i).CLDA(j).FixedBlocks)
                counter_CLDA=counter_CLDA+1;
                Perf_Plot(i).SuccessRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).SuccessRate;
                Perf_Plot(i).TravelTime(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).TravelTime;
                Perf_Plot(i).BitRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).FittsRate;
                Perf_Plot(i).Clicker(counter_CLDA)={PerfsDays(i).CLDA(j).FixedBlocks(k).Clicker};
                
            end
            
        end
        
    end
    
end

% figure for marking the days for sessions_based plots
figure('position',[400 400 1800 100]);
box on
set(gca,'FontSize',20)
vline(DayTrials)
yticks('')
xlim([0 length(Angle_x_AcrossTrials)])
xticks(DayTrials)
xtickangle(45)
xticklabels(Perf_Days)
%HighQualityFigsSVG('Fig_R2_2')

% figure for bit rates
figure('position',[400 400 1000 100]);
for i=1:length(Perf_Plot)
    
    plot(DayTrials(i),max(Perf_Plot(i).BitRate),'ok','MarkerFaceColor','k','MarkerSize',6)
    [m n]=max(Perf_Plot(i).BitRate);
    Index(i)=n;
    box off
    hold on  
end
yticks([0:0.2:0.8])
ylim([0,.8])
xlim([0 length(Angle_x_AcrossTrials)])
set(gca,'FontSize',10)
hold on
vline(DayTrials)
%HighQualityFigsSVG('Fig_R2_3')

[h,p,ci,stats]=ttest2([max(Perf_Plot(1).BitRate),max(Perf_Plot(2).BitRate),...
    max(Perf_Plot(3).BitRate)],...
[max(Perf_Plot(4).BitRate),max(Perf_Plot(5).BitRate)...
max(Perf_Plot(6).BitRate),max(Perf_Plot(7).BitRate),...
max(Perf_Plot(8).BitRate)])

effect_size=mean([max(Perf_Plot(1).BitRate),max(Perf_Plot(2).BitRate),...
    max(Perf_Plot(3).BitRate)])-mean...
([max(Perf_Plot(4).BitRate),max(Perf_Plot(5).BitRate)...
max(Perf_Plot(6).BitRate),max(Perf_Plot(7).BitRate),...
max(Perf_Plot(8).BitRate)])

% figure for success rate and travel time
fig=figure('position',[400 400 1000 100]);
left_color = [1 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
for i=1:length(Perf_Plot)
    
    plot(DayTrials(i),Perf_Plot(i).SuccessRate(Index(i)),'or','MarkerFaceColor','r','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,101])
xlim([0 length(Angle_x_AcrossTrials)])
xticks('')
set(gca,'FontSize',20)

yyaxis right
for i=1:length(Perf_Plot)
    
    plot(DayTrials(i),Perf_Plot(i).TravelTime(Index(i)),'ob','MarkerFaceColor','b','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,15])
set(gca,'FontSize',20)
xticks('')
hold on
vline(DayTrials)
%HighQualityFigsSVG('Fig_R2_4')

%%  Section 1-5: track the performance of PnP block by block for the based sessions chosen in section 1-4

clear all
close all
clc
%Experiment info: best sessions of PnP on each day
expts = [];
expts(end+1).yymmdd = '20200115';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'111243'};

expts(end+1).yymmdd = '20200117';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'113322'};

expts(end+1).yymmdd = '20200124';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'114404'};

expts(end+1).yymmdd = '20200127';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'142707'};

expts(end+1).yymmdd = '20200131';
expts(end).CLDA_hhmmss = {};
expts(end).perf_CenterOut_hhmmss = {};
expts(end).perf_RadialTask_hhmmss = {'111305'};

% Calculate the fitts rates, performance (success rates, time to target) 
%corresponding to blocks of PnP best sessions 
Day_Ref='20190521';

% for only PnP
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    PerfsDays(i).TestDate=yymmdd;
    
    NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
        [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);
    PerfsDays(i).TestDayCounter=NumDays;
    
    % going through the fixed blocks
    
    PerfsDays(i).PnP(1).SessionName='Plug&Play';
    FixedType='RadialTask';
    PerfsDays(i).FixedType=FixedType;
    FixedBlockCounter=0;
    for session=1:length(expt.perf_RadialTask_hhmmss)
        
        % go through datafiles in fixed sessions blocks
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer',yymmdd,FixedType,expt.perf_RadialTask_hhmmss{1,session},'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        
        % no matter of length to be 40, 80,...
        if mod(length(files),8)
            block_size=8;
            Block_Number=length(files)/8;
            
        elseif mod(length(files),10)
            block_size=10;
            Block_Number=length(files)/10;
        end
        
        load(fullfile(datadir,files(1).name));
        if isfield(TrialData,'ClickerState')
            Clicker=1;%'on';
            
            acc=0;fail=0;t=[];fail_index=[];correct_index=[];
            for Block=1:Block_Number
                FixedBlockCounter=FixedBlockCounter+1;
                
                for ii=(Block-1)*block_size+1:Block*block_size
                    
                    load(fullfile(datadir,files(ii).name));
                    if (length(TrialData.Time)-1)*1/TrialData.Params.UpdateRate < 20
                        if TrialData.TargetID ==  TrialData.SelectedTargetID
                            acc=acc+1;
                            correct_index = [correct_index ii];
                        elseif TrialData.TargetID ~=  TrialData.SelectedTargetID ...
                                && ~TrialData.SelectedTargetID == 0
                            fail = fail+1;
                            fail_index = [fail_index ii];
                        end
                        %t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
                        t=[t (length(TrialData.Time)-1)*1/TrialData.Params.UpdateRate];
                    end
                end
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).FittsRate=...
                    (3*max(0,acc-fail))/sum(t);
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).SuccessRate=...
                    acc/length(files)*100;
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).TravelTime=...
                    nanmean(t);
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).Clicker=...
                    Clicker;
                
            end
            
        elseif ~isfield(TrialData,'ClickerState')
            Clicker='off';
            
            acc=0;t=[];
            for Block=1:Block_Number
                FixedBlockCounter=FixedBlockCounter+1;
                
                for ii=(Block-1)*block_size+1:Block*block_size
                    load(fullfile(datadir,files(ii).name));
                    
                    if TrialData.ErrorID==0
                        acc=acc+1;
                    end
                    
                    t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
                end
                
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).FittsRate=...
                    (3*max(0,acc-fail))/sum(t);
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).SuccessRate=...
                    acc/length(files)*100;
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).TravelTime=...
                    nanmean(t);
                PerfsDays(i).PnP(1).FixedBlocks(FixedBlockCounter).Clicker=...
                    Clicker;              
            end       
        end    
    end   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the results
Perf_Days={};
for i=1:length(PerfsDays)
    
    Perf_Days(i)={num2str(PerfsDays(i).TestDayCounter)};
    counter=0;
    
    % it is PnP
    
    for k=1:length(PerfsDays(i).PnP(1).FixedBlocks)
        counter=counter+1;
        Perf_Plot(i).SuccessRate(counter)=PerfsDays(i).PnP(1).FixedBlocks(k).SuccessRate;
        Perf_Plot(i).TravelTime(counter)=PerfsDays(i).PnP(1).FixedBlocks(k).TravelTime;
        Perf_Plot(i).BitRate(counter)=PerfsDays(i).PnP(1).FixedBlocks(k).FittsRate;
        Perf_Plot(i).Clicker(counter)={PerfsDays(i).PnP(1).FixedBlocks(k).Clicker};
        
    end
     
end

% figure for bit rates
figure('position',[400 400 1000 300]);
Marker={'b','r','k','g','m'};
Points=1;
for i=1:length(Perf_Plot)
    plot(Points:Points-1+length(Perf_Plot(i).BitRate),...
        Perf_Plot(i).BitRate,'o-k','MarkerFaceColor',Marker{i},'MarkerSize',6)
    Points=Points+length(Perf_Plot(i).BitRate);
    %box off
    hold on
    
end

yticks([0:0.2:1])
ylim([0,1])
ylabel('Bits/ Sec')
xticks([1 3 5 8 10])
xticklabels({'day1','day2','day3','day4','day5'})
title('Block by block perfoamance for best sessions of PnP days')
set(gca,'FontSize',14)
%HighQualityFigs('Blocks_perf_PnP')

%% for section 2 start from begining
clear all
close all
clc


%% experiment info for long-term CLDA for section 2; three features: delta, beta, high gamma
expts = [];

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

expts(end+1).yymmdd = '20191122';
expts(end).CLDA_hhmmss = {'112414','141701'};
expts(end).perf_CenterOut_hhmmss = {'112414','141701','142524','143443'};
expts(end).perf_RadialTask_hhmmss = {'113743','114003','114800'};

% from this day using the new ModularBCI codes desgined by Daniel
expts(end+1).yymmdd = '20191127';
expts(end).CLDA_hhmmss = {'142835'};
expts(end).perf_CenterOut_hhmmss = {'142835'};
expts(end).perf_RadialTask_hhmmss = {'144302','144516','144837','145055','145243','145903'};

% RadialTypingMultiClick for the following day
expts(end+1).yymmdd = '20191203';
expts(end).CLDA_hhmmss = {'134453','135420'};
expts(end).perf_CenterOut_hhmmss = {'134453','135420'};
expts(end).perf_RadialTask_hhmmss = {'141002'};

% RadialTypingMultiClick for the following day
expts(end+1).yymmdd = '20191206';
expts(end).CLDA_hhmmss = {'110809'};
% CLDA blocks without fixed blocks: '140949','142301','142643'
expts(end).perf_CenterOut_hhmmss = {'110809'};
expts(end).perf_RadialTask_hhmmss = {'112558','113508','114445'};

% RadialTyping for the following day
expts(end+1).yymmdd = '20191217';
expts(end).CLDA_hhmmss = {'134525'};
expts(end).perf_CenterOut_hhmmss = {'134525'};
expts(end).perf_RadialTask_hhmmss = {'145444','145838','150222'};


%% Section 2-1: Reset experiment info after previous long-term CLDA (for section 2)

expts(end+1).yymmdd = '20200122';
expts(end).CLDA_hhmmss = {'113333','145510','151403'};
expts(end).perf_CenterOut_hhmmss = {'110924'};
expts(end).perf_RadialTask_hhmmss = {'152159'};

expts(end+1).yymmdd = '20200124';
expts(end).CLDA_hhmmss = {'140515'};
expts(end).perf_CenterOut_hhmmss = {'134519','135959','140515'};
expts(end).perf_RadialTask_hhmmss = {'141253'};

expts(end+1).yymmdd = '20200127';
expts(end).CLDA_hhmmss = {'114408','140423'};
expts(end).perf_CenterOut_hhmmss = {'112639','114408','140423'};
expts(end).perf_RadialTask_hhmmss = {'141809'};

expts(end+1).yymmdd = '20200131';
expts(end).CLDA_hhmmss = {'113822','133056'};
expts(end).perf_CenterOut_hhmmss = {'112800','134554'};
expts(end).perf_RadialTask_hhmmss = {'135215'};

expts(end+1).yymmdd = '20200204';
expts(end).CLDA_hhmmss = {'114712','133740'};
expts(end).perf_CenterOut_hhmmss = {'113550','133740'};
expts(end).perf_RadialTask_hhmmss = {'140041'};

%% Section 2-2: calculate and plot the angles for single Micro-CLDA and Reset KF and also the the average of 5 of them
close all
Sessions={};
% Considering one Ref for all calculation
Expt_Ref=expts(7); % the last day of CLDA of PnP
datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',Expt_Ref.yymmdd,'GangulyServer',...
    Expt_Ref.yymmdd,'CenterOut',Expt_Ref.CLDA_hhmmss{1,length(Expt_Ref.CLDA_hhmmss)},'BCI_CLDA');
files = dir(fullfile(datadir,'Data*.mat'));
load(fullfile(datadir,files(end).name));
C_Ref=TrialData.KalmanFilter{1,1}.C;

% calculate the angles relative to ref
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    Counter=0;
    for j=1:length(expt.CLDA_hhmmss)
        Sessions(end+1)={[yymmdd,'-',expt.CLDA_hhmmss{1,j}]};
        fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        if i<4 % following the previous structure of data recording
            datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
                'GangulyServer','Center-Out',yymmdd,expt.CLDA_hhmmss{1,j},'BCI_CLDA');
            
        else
            datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,'GangulyServer',...
                yymmdd,'CenterOut',expt.CLDA_hhmmss{1,j},'BCI_CLDA');
            
        end
        
        files = dir(fullfile(datadir,'Data*.mat'));
     
        for k=1:length(files)
            Counter=Counter+1;
            load(fullfile(datadir,files(k).name));
            C_New=TrialData.KalmanFilter{1,1}.C; 
            Angles_Plot(i).Xdirection(Counter)=subspace(C_New(:,3),C_Ref(:,3))*180/pi;
            Angles_Plot(i).Ydirection(Counter)=subspace(C_New(:,4),C_Ref(:,4))*180/pi;
        end
        
    end
    
        
end

% plot an example of one sample Micro-CLDA and One Reset KF
figure('position',[400 400 500 300]);

plot([Angles_Plot(1).Xdirection,Angles_Plot(11).Xdirection],'linewidth',2)
hold on
plot([Angles_Plot(1).Ydirection,Angles_Plot(11).Ydirection],'linewidth',2)
ylim([-1 90])
xlim([0 length([Angles_Plot(1).Xdirection,Angles_Plot(11).Xdirection])])
yticks(0:45:90)
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
set(gca,'FontSize',20)
xlabel('Trials')
box off
%HighQualityFigsSVG('Fig_R3_1')


% plot the average of 5 samples for Micro-CLDA and Reset KF
% for Micro-CLDA
MicroX=nan(7,48);
MicroY=nan(7,48);
for i=1:7
    if length(Angles_Plot(i).Xdirection)<48
        
        for j=1:length(Angles_Plot(i).Xdirection)
            MicroX(i,(48-1-length(Angles_Plot(i).Xdirection))+j)=Angles_Plot(i).Xdirection(j);
            MicroY(i,(48-1-length(Angles_Plot(i).Xdirection))+j)=Angles_Plot(i).Ydirection(j);  
        end
    else
        for j=1:length(Angles_Plot(i).Xdirection)
            MicroX(i,j)=Angles_Plot(i).Xdirection(j);
            MicroY(i,j)=Angles_Plot(i).Ydirection(j);
        end
    end
    
end

% for Reset KF
ResetX=nan(5,48);
ResetY=nan(5,48);
for i=8:12
    for j=1:length(Angles_Plot(i).Xdirection)
        ResetX(i-7,j)=Angles_Plot(i).Xdirection(j);
        ResetY(i-7,j)=Angles_Plot(i).Ydirection(j);   
    end
       
end

figure('position',[400 400 500 300]);
plot([nanmean(MicroX(:,1:48),1),nanmean(ResetX,1)],'linewidth',2)
hold on
plot([nanmean(MicroY(:,1:48),1),nanmean(ResetY,1)],'linewidth',2)
ylim([-1 90])
xlim([0 2*48])
yticks(0:45:90)
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
set(gca,'FontSize',20)
xlabel('Trials')
box off

%HighQualityFigsSVG('Fig_R3_2')

%% Section 2-3: plot KF convergence performance of Long-term CLDA in Center-Out with continuation of multiple Reset KF tests for Center-Out and RadialTask
% Calculate the fitts rates, performance (success rates, time to target) corresponding to CLDA blocks and also Reset KF tests

Day_Ref='20190521';

% for long-term CLDA & Reset KF
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    PerfsDays(i).TestDate=yymmdd;
    NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
        [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);
    PerfsDays(i).TestDayCounter=NumDays;
    

    
    if i<8 % for long-term CLDA days
        % for each day we only have one type of fixed block; so:
        FixedType='Center-Out';
        PerfsDays(i).FixedType=FixedType;
        
        % go througth the number of CLDA sessions
        for j=1:length(expt.CLDA_hhmmss)
            fprintf('CLDA:')
            fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
            PerfsDays(i).CLDA(j).SessionNameCLDA=expt.CLDA_hhmmss{1,j};
            
            if j<length(expt.CLDA_hhmmss)
                % go througth the number of related fixed sessions
                FixedBlockCounter=0;
                for jj=1:length(expt.perf_CenterOut_hhmmss)
                    
                    if str2num(expt.perf_CenterOut_hhmmss{1,jj})>= str2num(expt.CLDA_hhmmss{1,j}) && ...
                            str2num(expt.perf_CenterOut_hhmmss{1,jj})<str2num(expt.CLDA_hhmmss{1,j+1})
                        
                        FixedBlockCounter=FixedBlockCounter+1;
                        %calculate the type of fixed block, fitts rate & performance & clicker
                        %on/off
                        [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                        % save into structure:
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
                        
                    end
                    
                end
                
            elseif j==length(expt.CLDA_hhmmss)
                
                FixedBlockCounter=0;
                for jj=1:length(expt.perf_CenterOut_hhmmss)
                    
                    if str2num(expt.perf_CenterOut_hhmmss{1,jj})>= str2num(expt.CLDA_hhmmss{1,j})
                        
                        FixedBlockCounter=FixedBlockCounter+1;
                        %calculate the type of fixed block, fitts rate & performance & clicker
                        %on/off
                        [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                        % save into structure:
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                        PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
                        
                    end
                    
                end
                
            end
            
        end
        
    else % for Reset KF sessions
        
        % for each day we only have one type of fixed block; so:
        PerfsDays(i).FixedType='Combined:Center-Out&RadialTask';
        
        % go througth the number of CLDA sessions
        for j=[1,length(expt.CLDA_hhmmss)]
            k1=1; % for imagined blocks performance (before CLDA)
            k2=1; % for CLDA performance (after all CLDAs)
            fprintf('CLDA:')
            fprintf('\n%s-%s\n',yymmdd,expt.CLDA_hhmmss{1,j})
            
            % go througth the number of related fixed sessions
            FixedBlockCounter1=0;
            FixedBlockCounter2=0;
            
            % for centerout task
            for jj=1:length(expt.perf_CenterOut_hhmmss)
                FixedType='Center-Out';
                
                if str2num(expt.perf_CenterOut_hhmmss{1,jj})<str2num(expt.CLDA_hhmmss{1,1})
                    
                    PerfsDays(i).Imagined(k1).SessionName='Imagined_Perf';
                    FixedBlockCounter1=FixedBlockCounter1+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                    % save into structure:
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).FittsRate=FittsRate;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).SuccessRate=SuccessRate;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).TravelTime=TravelTime;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).Clicker=Clicker;
                    
                elseif str2num(expt.perf_CenterOut_hhmmss{1,jj})>=str2num(expt.CLDA_hhmmss{1,length(expt.CLDA_hhmmss)})
                    
                    PerfsDays(i).CLDA(k2).SessionNameCLDA=expt.CLDA_hhmmss{1,j};
                    FixedBlockCounter2=FixedBlockCounter2+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_CenterOut_hhmmss{1,jj},FixedType);
                    % save into structure:
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).FittsRate=FittsRate;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).SuccessRate=SuccessRate;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).TravelTime=TravelTime;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).Clicker=Clicker;
                    
                end
                
            end
            
            % for centerout task
            for jj=1:length(expt.perf_RadialTask_hhmmss)
                FixedType='RadialTask';
                
                if str2num(expt.perf_RadialTask_hhmmss{1,jj})<str2num(expt.CLDA_hhmmss{1,1})
                    
                    PerfsDays(i).Imagined(k1).SessionName='Imagined_Perf';
                    FixedBlockCounter1=FixedBlockCounter1+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_RadialTask_hhmmss{1,jj},FixedType);
                    % save into structure:
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).FittsRate=FittsRate;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).SuccessRate=SuccessRate;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).TravelTime=TravelTime;
                    PerfsDays(i).Imagined(k1).FixedBlocks(FixedBlockCounter1).Clicker=Clicker;
                    
                elseif str2num(expt.perf_RadialTask_hhmmss{1,jj})>=str2num(expt.CLDA_hhmmss{1,length(expt.CLDA_hhmmss)})
                    
                    PerfsDays(i).CLDA(k2).SessionNameCLDA=expt.CLDA_hhmmss{1,j};
                    FixedBlockCounter2=FixedBlockCounter2+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_RadialTask_hhmmss{1,jj},FixedType);
                    % save into structure:
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).FittsRate=FittsRate;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).SuccessRate=SuccessRate;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).TravelTime=TravelTime;
                    PerfsDays(i).CLDA(k2).FixedBlocks(FixedBlockCounter2).Clicker=Clicker;
                    
                end
                
            end
            
            
        end
        
    end
    
end

% plotting the results
Perf_Days={};

for i=1:length(PerfsDays)
    
    Perf_Days(i)={num2str(PerfsDays(i).TestDayCounter)};
    counter_CLDA=0;
    counter_Imagined=0;
    
    if length(PerfsDays(i).CLDA)==1
        if length(PerfsDays(i).CLDA.FixedBlocks)==1
            counter_CLDA=counter_CLDA+1;
            Perf_CLDA(i).SuccessRate(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks.SuccessRate;
            Perf_CLDA(i).TravelTime(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks.TravelTime;
            Perf_CLDA(i).BitRate(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks.FittsRate;
            Perf_CLDA(i).Clicker(counter_CLDA)={PerfsDays(i).CLDA.FixedBlocks.Clicker};
        else
            for k=1:length(PerfsDays(i).CLDA.FixedBlocks)
                counter_CLDA=counter_CLDA+1;
                Perf_CLDA(i).SuccessRate(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks(k).SuccessRate;
                Perf_CLDA(i).TravelTime(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks(k).TravelTime;
                Perf_CLDA(i).BitRate(counter_CLDA)=PerfsDays(i).CLDA.FixedBlocks(k).FittsRate;
                Perf_CLDA(i).Clicker(counter_CLDA)={PerfsDays(i).CLDA.FixedBlocks(k).Clicker};
            end   
        end
        
    else
        
        for j=1:length(PerfsDays(i).CLDA)
            
            if length(PerfsDays(i).CLDA(j).FixedBlocks)==1
                counter_CLDA=counter_CLDA+1;
                Perf_CLDA(i).SuccessRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks.SuccessRate;
                Perf_CLDA(i).TravelTime(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks.TravelTime;
                Perf_CLDA(i).BitRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks.FittsRate;
                Perf_CLDA(i).Clicker(counter_CLDA)={PerfsDays(i).CLDA(j).FixedBlocks.Clicker};
                
            else
                for k=1:length(PerfsDays(i).CLDA(j).FixedBlocks)
                    counter_CLDA=counter_CLDA+1;
                    Perf_CLDA(i).SuccessRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).SuccessRate;
                    Perf_CLDA(i).TravelTime(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).TravelTime;
                    Perf_CLDA(i).BitRate(counter_CLDA)=PerfsDays(i).CLDA(j).FixedBlocks(k).FittsRate;
                    Perf_CLDA(i).Clicker(counter_CLDA)={PerfsDays(i).CLDA(j).FixedBlocks(k).Clicker};
                end
                   
            end   
        end
        
    end
    
    if ~isempty(PerfsDays(i).Imagined)
        
        if length(PerfsDays(i).Imagined.FixedBlocks)==1
            counter_Imagined=counter_Imagined+1;
            Perf_Imagined(i).SuccessRate(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks.SuccessRate;
            Perf_Imagined(i).TravelTime(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks.TravelTime;
            Perf_Imagined(i).BitRate(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks.FittsRate;
            Perf_Imagined(i).Clicker(counter_Imagined)={PerfsDays(i).Imagined.FixedBlocks.Clicker};
        else
            for k=1:length(PerfsDays(i).Imagined.FixedBlocks)
                counter_Imagined=counter_Imagined+1;
                Perf_Imagined(i).SuccessRate(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks(k).SuccessRate;
                Perf_Imagined(i).TravelTime(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks(k).TravelTime;
                Perf_Imagined(i).BitRate(counter_Imagined)=PerfsDays(i).Imagined.FixedBlocks(k).FittsRate;
                Perf_Imagined(i).Clicker(counter_Imagined)={PerfsDays(i).Imagined.FixedBlocks(k).Clicker};
            end
        end  
    end
    
end

% plot the corresponding bitrate for one single sample
figure('position',[400 400 400 250]);
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#21409A';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
bar(1,mean(Perf_CLDA(1).BitRate),0.3,'FaceColor',color,'EdgeColor',color)
hold on

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#FFDE17';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
bar([2,3],[mean(Perf_Imagined(11).BitRate),...
    mean(Perf_CLDA(11).BitRate)],0.3,'FaceColor',color,'EdgeColor',color)

yticks([0.2 .4 .6])
xticks([1 2 3])
xlim([0.5 3.5])
set(gca,'FontSize',20)
box off
%HighQualityFigsSVG('Fig_R3_3')

% plot the corresponding bitrate for all samples
figure('position',[400 400 400 250]);
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#21409A';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
bar(1,mean([Perf_CLDA(1).BitRate,...
    Perf_CLDA(2).BitRate,...
    Perf_CLDA(3).BitRate,...
    Perf_CLDA(4).BitRate,...
    Perf_CLDA(5).BitRate,...
    Perf_CLDA(6).BitRate,...
    Perf_CLDA(7).BitRate]),0.3,'FaceColor',color,'EdgeColor',color)
hold on
errorbar(1,mean([Perf_CLDA(1).BitRate,...
    Perf_CLDA(2).BitRate,...
    Perf_CLDA(3).BitRate,...
    Perf_CLDA(4).BitRate,...
    Perf_CLDA(5).BitRate,...
    Perf_CLDA(6).BitRate,...
    Perf_CLDA(7).BitRate]),std([Perf_CLDA(1).BitRate,...
    Perf_CLDA(2).BitRate,...
    Perf_CLDA(3).BitRate,...
    Perf_CLDA(4).BitRate,...
    Perf_CLDA(5).BitRate,...
    Perf_CLDA(6).BitRate,...
    Perf_CLDA(7).BitRate]))
hold on
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#FFDE17';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
bar([2 3],[mean([Perf_Imagined(8).BitRate,...
    Perf_Imagined(9).BitRate,...
    Perf_Imagined(10).BitRate,...
    Perf_Imagined(11).BitRate,...
    Perf_Imagined(12).BitRate]),...
    mean([Perf_CLDA(8).BitRate,...
    Perf_CLDA(9).BitRate,...
    Perf_CLDA(10).BitRate,...
    Perf_CLDA(11).BitRate,...
    Perf_CLDA(12).BitRate])],0.3,'FaceColor',color,'EdgeColor',color)
hold on
errorbar([2 3],[mean([Perf_Imagined(8).BitRate,...
    Perf_Imagined(9).BitRate,...
    Perf_Imagined(10).BitRate,...
    Perf_Imagined(11).BitRate,...
    Perf_Imagined(12).BitRate]),...
    mean([Perf_CLDA(8).BitRate,...
    Perf_CLDA(9).BitRate,...
    Perf_CLDA(10).BitRate,...
    Perf_CLDA(11).BitRate,...
    Perf_CLDA(12).BitRate])],[std([Perf_Imagined(8).BitRate,...
    Perf_Imagined(9).BitRate,...
    Perf_Imagined(10).BitRate,...
    Perf_Imagined(11).BitRate,...
    Perf_Imagined(12).BitRate]),...
    std([Perf_CLDA(8).BitRate,...
    Perf_CLDA(9).BitRate,...
    Perf_CLDA(10).BitRate,...
    Perf_CLDA(11).BitRate,...
    Perf_CLDA(12).BitRate])])

yticks([0.2 .4 .6 0.8])
ylim([0 0.9])
xticks([1 2 3])
xlim([0.5 3.5])
set(gca,'FontSize',20)
box off
%HighQualityFigsSVG('Fig_R3_4')


%% Section 2-4: plot the trend of changes across CLDA trials for feature weights after Reset KF 
% plot across array
ch_layout = [
    96	84	76	95	70	82	77	87	74	93	66	89	86	94	91	79
    92	65	85	83	68	75	78	81	72	69	88	71	80	73	90	67
    62	37	56	48	43	44	60	33	49	64	58	59	63	61	51	34
    45	53	55	52	35	57	38	50	54	39	47	42	36	40	46	41
    19	2	10	21	30	23	17	28	18	1	8	15	32	27	9	3
    24	13	6	4	7	16	22	5	20	14	11	12	29	26	31	25
    124	126	128	119	110	113	111	122	117	125	112	98	104	116	103	106
    102	109	99	101	121	127	105	120	107	123	118	114	108	115	100	97];
[R,Col] = size(ch_layout);
Nch = 128;
feature_strs={'delta','beta','hg'};

% running examples to observe the plots and choose the trials for next figures
i=11;
expt = expts(i);
yymmdd = expt.yymmdd;

for j=1:length(expt.CLDA_hhmmss)
    
    % go through datafiles in CLDA blocks
    if i<4 % following the previous structure of data recording
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
    else
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,'GangulyServer',...
            yymmdd,'CenterOut',expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
    end
          
    files = dir(fullfile(datadir,'Data*.mat'));
    
    for kk=1:ceil(length(files)/16)
        figure('position',[1200 0 400 1300]);
        suptitle(['CLDA',num2str(j),'-',num2str(kk),'; after reset-KF'])
        
        for k=(1+16*(kk-1)):min(16*kk,length(files))
            load(fullfile(datadir,files(k).name));
            C_New=TrialData.KalmanFilter{1,1}.C;
            Cx = C_New(:,3);
            Cy = -C_New(:,4);
            cc = brewermap(75,'YlOrRd');
            
            for feature=1:3
                Cx_feature = Cx(128*(feature-1)+1:128*feature,1);
                Cy_feature = Cy(128*(feature-1)+1:128*feature,1);
                Cvec=[Cx_feature,Cy_feature];
                for ii=1:size(Cvec,1)
                    mag(ii) = norm(Cvec(ii,:));
                end
                
                % put into plotting matrix
                x = zeros(size(ch_layout));
                y = zeros(size(ch_layout));
                M = zeros(size(ch_layout));
                for ch=1:Nch
                    [r,c] = find(ch_layout == ch);
                    x(r,c) = c;
                    y(r,c) = (R-r)+1;
                    M(r,c) = mag(ch);
                end
                
                % plot
                if mod(k,16)==0
                    ax=subplot(16,3,3*(16-1)+feature);
                else
                    ax=subplot(16,3,3*(mod(k,16)-1)+feature);
                end
                
                imagesc(M)
                colormap(cc)
                set(ax,'XTick',[],'YTick',[])
                %colorbar(ax,'location','westoutside')
                title(sprintf('%s',feature_strs{feature}))
                
            end
            
        end
       
    end
    
end

% running the selected trials for figures
% the trials chosed for Reset KF from day 11: 1,11,35,48: so:
% j=1 k=1,11 
% j=2 K=10,24
% the trials chosed for Micro-CLDA from day 1: 1,11,35,48; so:
% j=1 k=1,11 
% j=3 K=3,16

% each trial will be run seperately to generate the figures
figure('position',[400 400 450 80]);
i=1;
expt = expts(i);
yymmdd = expt.yymmdd;
for j=3
    
    % go through datafiles in CLDA blocks
    if i<4 % following the previous structure of data recording
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
    else
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,'GangulyServer',...
            yymmdd,'CenterOut',expt.CLDA_hhmmss{1,j},'BCI_CLDA');
        
    end
    
    files = dir(fullfile(datadir,'Data*.mat'));
    
    for k=16
        load(fullfile(datadir,files(k).name));
        C_New=TrialData.KalmanFilter{1,1}.C;
        Cx = C_New(:,3);
        Cy = -C_New(:,4);
        cc = brewermap(75,'YlOrRd');
        
        for feature=1:3
            Cx_feature = Cx(128*(feature-1)+1:128*feature,1);
            Cy_feature = Cy(128*(feature-1)+1:128*feature,1);
            Cvec=[Cx_feature,Cy_feature];
            for ii=1:size(Cvec,1)
                mag(ii) = norm(Cvec(ii,:));
            end
            
            % put into plotting matrix
            x = zeros(size(ch_layout));
            y = zeros(size(ch_layout));
            M = zeros(size(ch_layout));
            for ch=1:Nch
                [r,c] = find(ch_layout == ch);
                x(r,c) = c;
                y(r,c) = (R-r)+1;
                M(r,c) = mag(ch);
            end
            
            % plot
            ax=subplot(1,3,feature);
            imagesc(M)
            colormap(cc)
            set(ax,'XTick',[],'YTick',[])
            %colorbar(ax,'location','westoutside')
            title(sprintf('%s',feature_strs{feature}))
            
        end
    end
end
