%% 
% Reza note: The calculated Fitt's rate for the last day (5 points):
% This is the value for the last day of ltCLDA for comparison: Mean=0.3958 +-(0.0489)

%% plug_and_play sessions
clear, clc, close all

%% experiment info for only Decoder (A) which was Original stable decoder
expts = [];

expts(end+1).yymmdd = '20190724';
expts(end).perf_hhmmss = {'111324','134612'};

expts(end+1).yymmdd = '20190725';
expts(end).perf_hhmmss = {'111205','134629'};

expts(end+1).yymmdd = '20190726';
expts(end).perf_hhmmss = {'111631'};

expts(end+1).yymmdd = '20190730';
expts(end).perf_hhmmss = {'111316'};


%% My own code: Calculate the fitts rates, performance (success rates, time to target)

Day_Ref='20190626';

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    PerfsDays(i).TestDate=yymmdd;
    
    FixedType='Center-Out';
    PerfsDays(i).FixedType=FixedType;
    
    NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
        [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);
    PerfsDays(i).TestDayCounter=NumDays;
    
    % go througth the number of related fixed sessions
    FixedBlockCounter=0;
    for jj=1:length(expt.perf_hhmmss)
        FixedBlockCounter=FixedBlockCounter+1;
        %calculate the type of fixed block, fitts rate & performance & clicker
        %on/off
        [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_hhmmss{1,jj},FixedType);
        
        % save into structure:
        PerfsDays(i).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
        PerfsDays(i).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
        PerfsDays(i).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
        PerfsDays(i).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
           
    end
      
end

% plotting the results
Perf_Days={};

for i=1:length(PerfsDays)
    
    Perf_Days(i)={num2str(PerfsDays(i).TestDayCounter)};
    counter=0;
    
    for k=1:length(PerfsDays(i).FixedBlocks)
        counter=counter+1;
        Perf_Plot(i).SuccessRate(counter)=PerfsDays(i).FixedBlocks(k).SuccessRate;
        Perf_Plot(i).TravelTime(counter)=PerfsDays(i).FixedBlocks(k).TravelTime;
        Perf_Plot(i).BitRate(counter)=PerfsDays(i).FixedBlocks(k).FittsRate;
        Perf_Plot(i).Clicker(counter)={PerfsDays(i).FixedBlocks(k).Clicker};
        
    end
    
end

% figure for bit rates
figure('position',[445 326 400 110]);
for i=1:length(Perf_Plot)
    
    plot(i-0.5,mean(Perf_Plot(i).BitRate),'ok','MarkerFaceColor','k','MarkerSize',6)
    box off
    hold on  
end
yticks(0:0.2:0.6)
ylim([0 0.6])
xlim([0,4])
set(gca,'FontSize',12)
hold on
%HighQualityFigsSVG('Fig_R1_5')
% shadow area
x = [0 4 4 0];
y = [0.3958-(2*0.0489) 0.3958-(2*0.0489) 0.3958+(2*0.0489) 0.3958+(2*0.0489)];
patch(x,y,'b','FaceAlpha',.3,'EdgeAlpha',.01)


% figure for success rate and travel time
%fig=figure('position',[400 400 1800 300]);
fig=figure('position',[445 326 400 110]);
left_color = [1 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
for i=1:length(Perf_Plot)
    
    plot(i-0.5,mean(Perf_Plot(i).SuccessRate),'or','MarkerFaceColor','r','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,105])
xlim([0,4])
set(gca,'FontSize',12)

yyaxis right
for i=1:length(Perf_Plot)
    
    plot(i-0.5,mean(Perf_Plot(i).TravelTime),'ob','MarkerFaceColor','b','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,15])
set(gca,'FontSize',12)

hold on
%vline(SessionsNum)
%HighQualityFigsSVG('Fig_R1_6')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% collect beh measures per target

target_ang = 0:45:360-45;

SR      = cell(1,8); % Success Rate
TtST	= cell(1,8); % Time to Start Target
TtRT	= cell(1,8); % Time to Reach Target
RR      = cell(1);   % Reward Rate (SR / (TtST + TtRT))
MD      = cell(1,8); % Max Deviation from Optimal Trajectory
Pr      = cell(1,8); % Avg. Progress
PL      = cell(1,8); % Path Length
Jerk	= cell(1,8);
Perf    = {};

% go through expts
ct = 0;
success = 1;
success_prev = 1;
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    for ii=1:length(expt.hhmmss),
        hhmmss = expt.hhmmss{ii};
        ct = ct + 1;
        
        % per session
        SR_sess     = cell(1,8); % Success Rate
        TT_sess     = cell(1,8); % Total Targets
        TtST_sess   = cell(1,8); % Time to Start Target
        TtRT_sess   = cell(1,8); % Time to Reach Target
        MD_sess     = cell(1,8); % Max Deviation from Optimal Trajectory
        Pr_sess     = cell(1,8); % Avg. Progress
        PL_sess     = cell(1,8); % Path Length
        Jerk_sess   = cell(1,8); % Jerk
        
        % go through datafiles in fixed blocks
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        fprintf('\n%s-%s\n',yymmdd,hhmmss)
        files = dir(fullfile(datadir,'Data*.mat'));
        
        for iii=9:length(files), % skip first 8 trials
        
%         for iii=1:length(files),
            % load data
            load(fullfile(datadir,files(iii).name));
            target_idx = TrialData.TargetAngle==target_ang;
            
            % Compute Behavioral Metrics
            
            % Success Rate
            SR_sess{target_idx}(end+1) = TrialData.ErrorID==0;
            success_prev = success;
            success = TrialData.ErrorID==0;
            
            % Total Targets
            if TrialData.Params.CenterReset || success_prev==0, % no start target
                switch TrialData.ErrorID, % error during trial
                    case 0, TT_sess{target_idx}(end+1) = 1;
                    otherwise, TT_sess{target_idx}(end+1) = 0;
                end
            else, % include start target
                switch TrialData.ErrorID, % error during trial
                    case 0, TT_sess{target_idx}(end+1) = 2;
                    case 1, TT_sess{target_idx}(end+1) = 0;
                    otherwise, TT_sess{target_idx}(end+1) = 1;
                end
            end

            % Time to Start Target
            % skip going to start target if 1st trial of block or if no
            % start target pd.
            if TrialData.Params.CenterReset || ...
                    mod(TrialData.Trial,TrialData.Params.NumTrialsPerBlock)==1 || ...
                    success_prev==0,
                TtST_sess{target_idx}(end+1) = NaN;
            else,
                idx = find(strcmp({TrialData.Events.Str},'Start Target'));
                if TrialData.ErrorID==1, % error getting to start target
                    TtST_sess{target_idx}(end+1) = TrialData.Params.MaxStartTime;
                else,
                    TtST_sess{target_idx}(end+1) = TrialData.Events(idx+1).Time-...
                        TrialData.Events(idx).Time;
                end
            end
            
            % Time to Reach Target
            if TrialData.ErrorID==0,
                idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                TtRT_sess{target_idx}(end+1) = TrialData.Time(end)-...
                    TrialData.Events(idx).Time;
            elseif TrialData.ErrorID==1,
                TtRT_sess{target_idx}(end+1) = NaN;
            else,
                TtRT_sess{target_idx}(end+1) = TrialData.Params.MaxReachTime;
            end
            
            % Progress
            if 0
                if TrialData.ErrorID==0,
                    idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                    tidx = TrialData.Time >= TrialData.Events(idx).Time;
                    dist = point_to_line_distance(...
                        TrialData.CursorState(1:2,tidx)',...
                        TrialData.CursorState(1:2,1)',TrialData.TargetPosition);
                    MD{target_idx}(end+1) = max(dist);
                else,
                    MD{target_idx}(end+1) = NaN;
                end
            end
            % Max Deviation from Optimal Traj during Reach Target
            if TrialData.ErrorID==0,
                idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                tidx = TrialData.Time >= TrialData.Events(idx).Time;
                V = TrialData.CursorState(3:4,tidx)';
                Vopt = TrialData.IntendedCursorState(3:4,tidx)';
                proj_v_vopt = Vopt .* dot(V,Vopt,2) ./ dot(Vopt,Vopt,2);
                progress = sqrt(dot(proj_v_vopt,proj_v_vopt,2));
                
                Pr_sess{target_idx}(end+1) = mean(progress);
            else,
                Pr_sess{target_idx}(end+1) = NaN;
            end
            
            % Normalized Path Length
            if TrialData.ErrorID==0,
                idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                tidx = TrialData.Time >= TrialData.Events(end).Time;
                P = TrialData.CursorState(1:2,tidx)';
                dP = diff(P,1,1);
                PL_sess{target_idx}(end+1) = sum(sqrt(sum(P.^2))) / ...
                    TrialData.Params.ReachTargetRadius;
            else,
                PL_sess{target_idx}(end+1) = NaN;
            end
            
            % Jerk
            if TrialData.ErrorID==0,
                idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                tidx = TrialData.Time >= TrialData.Events(end).Time;
                xpos = TrialData.CursorState(1,tidx);
                ypos = TrialData.CursorState(2,tidx);
                xvel = TrialData.CursorState(3,tidx);
                yvel = TrialData.CursorState(4,tidx);
                xacc = gradient(xvel,TrialData.Time(tidx));
                yacc = gradient(yvel,TrialData.Time(tidx));
                xjrk = gradient(xacc,TrialData.Time(tidx));
                yjrk = gradient(yacc,TrialData.Time(tidx));
                Jerk_sess{target_idx}(end+1) = (mean(sqrt(xjrk.^2 + yjrk.^2)));
            else,
                Jerk_sess{target_idx}(end+1) = NaN;
            end
            
        end
        
        % Per Session Measures
        AvgRR_sess = 60*sum(cat(2,TT_sess{:})) / ...
            (nansum(cat(2,TtST_sess{:})) + nansum(cat(2,TtRT_sess{:})));
        AvgSR_sess = mean(cat(2,SR_sess{:}));
        AvgTtST = nanmean(cat(2,TtST_sess{:}));
        AvgTtRT = nanmean(cat(2,TtRT_sess{:}));
        ITR_RT_sess = ...
            (1/AvgTtRT) * ...
            log2((TrialData.Params.ReachTargetRadius ...
            + TrialData.Params.TargetSize)/TrialData.Params.TargetSize);
        ITR_sess = ...
            (1/(nanmean([AvgTtST,AvgTtRT]))) * ...
            log2((TrialData.Params.ReachTargetRadius ...
            + TrialData.Params.TargetSize)/TrialData.Params.TargetSize);
        AvgJerk_sess = nanmean(cat(2,Jerk_sess{:}));
        
        % Output
        fprintf('Trials: %i\n',length(files))
        fprintf('Target Radius: %i\n',TrialData.Params.ReachTargetRadius)
        fprintf('Center-Reset: %i\n',TrialData.Params.CenterReset)
        fprintf('A: %.03f\n',TrialData.Params.KF.A(3,3))
        fprintf('W: %i\n',TrialData.Params.KF.W(3,3))
        fprintf('G: %i\n',TrialData.Params.Gain)
        try
            fprintf('Mask: %s\n',expt.masks{ii})
        end
        fprintf('Success Rate: %i%%\n',round(100*AvgSR_sess))
        fprintf('Time to Target: %.1f sec\n',nanmean([AvgTtST,AvgTtRT]))
        fprintf('Fitt''s IT Rate: %.2f bits/sec\n',ITR_sess)
        fprintf('Smoothness: %.2f \n',AvgJerk_sess)
        
        % store measures in cell array
        Perf{ct,1} = sprintf('%s-%s',yymmdd,hhmmss);
        Perf{ct,2} = length(files);
        Perf{ct,3} = round(100*AvgSR_sess);
        Perf{ct,4} = round(AvgRR_sess,3);
        Perf{ct,5} = round(AvgTtST,3);
        Perf{ct,6} = round(AvgTtRT,3);
        Perf{ct,7} = round(nanmean([AvgTtST,AvgTtRT]),3);
        Perf{ct,8} = round(ITR_RT_sess,3);
        Perf{ct,9} = round(ITR_sess,3);
        Perf{ct,10} = round(AvgJerk_sess,3);
        
    end
end

% save Performance in csv
PerformanceTable = cell2table(Perf);
PerformanceTable.Properties.VariableNames = {
    'Session','Trials','SuccessRate',...
    'RewardRate','TimeToStartTarget','TimeToReachTarget','AvgTimeToTarget',...
    'InfoTransferRateReachOnly','InfoTransferRateAll','Jerk'};
writetable(PerformanceTable,'PlugAndPlayExperimentPerformance.csv')

%% timeseries plot: plug and play vs. session
expts = [];

expts(end+1).yymmdd = '20190724';
expts(end).hhmmss = {'111324','134612','140804'};
expts(end).controller = {1,1,2};
expts(end).elapsedTime = {0,1,0};

expts(end+1).yymmdd = '20190725';
expts(end).hhmmss = {'105857','111205','112929','134629','140548'};
expts(end).controller = {2,1,3,1,4};
expts(end).elapsedTime = {1,2,0,3,0};

expts(end+1).yymmdd = '20190726';
expts(end).hhmmss = {'103925','105253','110243','111631','113418'};
expts(end).controller = {4,3,2,1,5};
expts(end).elapsedTime = {1,2,3,4,0};

% expts(end+1).yymmdd = '20190730';
% expts(end).hhmmss = {'111316','112918'};
% expts(end).controller = {1,6};
% expts(end).elapsedTime = {10,0};


% collect beh measures per target

X       = cell(1,11); % Controller (for x-axis of plot)
SR      = cell(1,11); % Success Rate
TtST	= cell(1,11); % Time to Start Target
TtRT	= cell(1,11); % Time to Reach Target
ITR     = cell(1,11); % Fitt's ITR

% go through expts
ct = 0;
success = 1;
success_prev = 1;
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    for ii=1:length(expt.hhmmss),
        hhmmss = expt.hhmmss{ii};
        ct = ct + 1;
        
        % per session
        SR_sess     = []; % Success Rate
        TtST_sess   = []; % Time to Start Target
        TtRT_sess   = []; % Time to Reach Target
        
        % go through datafiles in fixed blocks
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        fprintf('\n%s-%s\n',yymmdd,hhmmss)
        files = dir(fullfile(datadir,'Data*.mat'));
        for iii=9:length(files), % skip first 8 trials
        
%         for iii=1:length(files),
            % load data
            load(fullfile(datadir,files(iii).name));
            
            % Compute Behavioral Metrics
            
            % Success Rate
            SR_sess(end+1) = TrialData.ErrorID==0;
            success_prev = success;
            success = TrialData.ErrorID==0;
            
            % Time to Start Target
            % skip going to start target if 1st trial of block or if no
            % start target pd.
            if TrialData.Params.CenterReset || ...
                    mod(TrialData.Trial,TrialData.Params.NumTrialsPerBlock)==1 || ...
                    success_prev==0,
                TtST_sess(end+1) = NaN;
            else,
                idx = find(strcmp({TrialData.Events.Str},'Start Target'));
                if TrialData.ErrorID==1, % error getting to start target
                    TtST_sess(end+1) = TrialData.Params.MaxStartTime;
                else,
                    TtST_sess(end+1) = TrialData.Events(idx+1).Time-...
                        TrialData.Events(idx).Time;
                end
            end
            
            % Time to Reach Target
            if TrialData.ErrorID==0,
                idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                TtRT_sess(end+1) = TrialData.Time(end)-...
                    TrialData.Events(idx).Time;
            elseif TrialData.ErrorID==1,
                TtRT_sess(end+1) = NaN;
            else,
                TtRT_sess(end+1) = TrialData.Params.MaxReachTime;
            end
            
        end
        
        % Per Session Measures
        idx = expt.elapsedTime{ii} + 1;
        X{idx}(end+1)    = idx-1;
        SR{idx}(end+1)   = mean(SR_sess);
        TtST{idx}(end+1) = nanmean(TtST_sess);
        TtRT{idx}(end+1) = nanmean(TtRT_sess);
        ITR{idx}(end+1)  = ...
            (1/(nanmean([TtST_sess,TtRT_sess]))) * ...
            log2((TrialData.Params.ReachTargetRadius ...
            + TrialData.Params.TargetSize)/TrialData.Params.TargetSize);
        
    end
end

% some stats - t-test 0 vs not 0 and regress 1:5
[h1,p1] = ttest2(cat(2,SR{1})',cat(2,SR{2:end})');
mu1 = mean(cat(2,SR{2:end})) - mean(cat(2,SR{1}));
[b1,~,~,~,S1] = regress(cat(2,SR{2:end})',[ones(8,1),cat(2,X{2:end})']);
fprintf('\nSR:\n')
fprintf('t-test - eff_sz=%.03f,   p=%.03f\n',mu1,p1)
fprintf('regression - slope=%.03f, r2=%.03f, p=%.03f\n',b1(2),S1(1),S1(3))

TtT0 = mean(cat(1,cat(2,TtST{1}),cat(2,TtRT{1})));
TtT1 = mean(cat(1,cat(2,TtST{2:end}),cat(2,TtRT{2:end})));
[h2,p2] = ttest2(cat(2,TtT0)',cat(2,TtT1)');
mu2 = mean(cat(2,TtT1)) - mean(cat(2,TtT0));
[b2,~,~,~,S2] = regress(TtT1',[ones(8,1),cat(2,X{2:end})']);
fprintf('\nTtT:\n')
fprintf('t-test - eff_sz=%.03f,   p=%.03f\n',mu2,p2)
fprintf('regression - slope=%.03f, r2=%.03f, p=%.03f\n',b2(2),S2(1),S2(3))

[h3,p3] = ttest2(cat(2,ITR{1})',cat(2,ITR{2:end})');
mu3 = mean(cat(2,ITR{2:end})) - mean(cat(2,ITR{1}));
[b3,bint3,~,~,S3] = regress(cat(2,ITR{2:end})',[ones(8,1),cat(2,X{2:end})']);
fprintf('\nITR:\n')
fprintf('t-test - eff_sz=%.03f,   p=%.03f\n',mu3,p3)
fprintf('regression - int=%.03f, slope=%.03f, CI: [%.03f, %.03f]\nr2=%.03f, p=%.03f\n',...
    b3(1),b3(2),bint3(2,1),bint3(2,2),S3(1),S3(3))

% scatter plot
figure('position',[672 695 519 276]);
cc = get(groot,'defaultAxesColorOrder');

subplot(2,1,1)

yyaxis left, hold on
scatter(cat(2,X{:})-.1,cat(2,SR{:}),'o','markeredgecolor','w',...
    'sizedata',30,'markerfacecolor',cc(1,:))
plot([1,4]-.1,b1(1)+b1(2)*[1,4],'-','color',cc(1,:))
ylabel('% Success')
title('BCI Performance')
xlim([-0.5,4.5])
ylim([0,1])

TtT = mean(cat(1,cat(2,TtST{:}),cat(2,TtRT{:})));
yyaxis right, hold on
scatter(cat(2,X{:})+.1,TtT,'o','markeredgecolor','w',...
    'sizedata',30,'markerfacecolor',cc(2,:))
plot([1,4]+.1,b2(1)+b2(2)*[1,4],'-','color',cc(2,:))
ylabel('Time to Target')
title('BCI Performance')
xlim([-0.5,4.5])
ylim([0,15])
set(gca,'XTick',0:4,'XTickLabel',[])

subplot(2,1,2), hold on
scatter(cat(2,X{:}),cat(2,ITR{:}),'o','markeredgecolor','w',...
    'sizedata',30,'markerfacecolor','k')
plot([1,4],b3(1)+b3(2)*[1,4],'-','color','k')
ylabel('bit / sec')
xlabel('Days Since Decoder was Trained')
title('Fitt''s ITR')
xlim([-0.5,4.5])
ylim([0.1,0.55])
set(gca,'XTick',0:4,'XTickLabel',0:.5:4,'YTick',[.2,.4])


% save fig
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/PlugAndPlay/PerformanceDegradatation.pdf',...
%     '-pdf','-png','-r300','-painters')
% close(gcf)


%% timeseries plot: plug and play vs. session

expts = [];

expts(end+1).yymmdd = '20190724';
expts(end).hhmmss = {'111324','134612','140804'};
expts(end).controller = {1,1,2};
% expts(end).elapsedTime = {0,1,0};
expts(end).elapsedTime = {1,2,2};

expts(end+1).yymmdd = '20190725';
expts(end).hhmmss = {'105857','111205','112929','134629','140548'};
expts(end).controller = {2,1,3,1,4};
% expts(end).elapsedTime = {1,2,0,3,0};
expts(end).elapsedTime = {3,3,3,4,4};

expts(end+1).yymmdd = '20190726';
expts(end).hhmmss = {'103925','105253','110243','111631','113418'};
expts(end).controller = {4,3,2,1,5};
% expts(end).elapsedTime = {1,2,3,4,0};
expts(end).elapsedTime = {5,5,5,5,5};



% collect beh measures per target
X       = cell(1,5); % Controller (for x-axis of plot)
SR      = cell(1,5); % Success Rate
TtST	= cell(1,5); % Time to Start Target
TtRT	= cell(1,5); % Time to Reach Target
ITR     = cell(1,5); % Fitt's ITR

% go through expts
for decoder=1:5,
    for i=1:length(expts),
        expt = expts(i);
        yymmdd = expt.yymmdd;
        for ii=1:length(expt.hhmmss),
            hhmmss = expt.hhmmss{ii};
            
            % skip if session is not using decoder
            if expt.controller{ii}~=decoder,
                continue;
            end

            % per session
            SR_sess     = []; % Success Rate
            TtST_sess   = []; % Time to Start Target
            TtRT_sess   = []; % Time to Reach Target

            % go through datafiles in fixed blocks
            datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
                'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
            fprintf('\n%s-%s\n',yymmdd,hhmmss)
            files = dir(fullfile(datadir,'Data*.mat'));
            
            success = 1;
            for iii=1:length(files),
                % load data
                load(fullfile(datadir,files(iii).name));

                % Compute Behavioral Metrics

                % Success Rate
                success_prev = success;
                success = TrialData.ErrorID==0;
                SR_sess(end+1) = success;
                
                % Time to Start Target
                % skip going to start target if 1st trial of block or if no
                % start target pd.
                if TrialData.Params.CenterReset || ...
                        mod(TrialData.Trial,TrialData.Params.NumTrialsPerBlock)==1 || ...
                        success_prev==0,
                    TtST_sess(end+1) = NaN;
                else,
                    idx = find(strcmp({TrialData.Events.Str},'Start Target'));
                    if TrialData.ErrorID==1, % error getting to start target
                        TtST_sess(end+1) = TrialData.Params.MaxStartTime;
                    else,
                        TtST_sess(end+1) = TrialData.Events(idx+1).Time-...
                            TrialData.Events(idx).Time;
                    end
                end

                % Time to Reach Target
                if TrialData.ErrorID==0,
                    idx = find(strcmp({TrialData.Events.Str},'Reach Target'));
                    TtRT_sess(end+1) = TrialData.Time(end)-...
                        TrialData.Events(idx).Time;
                elseif TrialData.ErrorID==1,
                    TtRT_sess(end+1) = NaN;
                else,
                    TtRT_sess(end+1) = TrialData.Params.MaxReachTime;
                end

            end

            % Per Session Measures
            idx                  = expt.elapsedTime{ii} + 1;
            X{decoder}(end+1)    = idx-1;
            SR{decoder}(end+1)   = mean(SR_sess);
            TtST{decoder}(end+1) = nanmean(TtST_sess);
            TtRT{decoder}(end+1) = nanmean(TtRT_sess);
            ITR{decoder}(end+1)  = ...
                (1/(nanmean([TtST_sess,TtRT_sess]))) * ...
                log2((TrialData.Params.ReachTargetRadius ...
                + TrialData.Params.TargetSize)/TrialData.Params.TargetSize);

        end % sessions
    end % days
end % decoder


% Plot for performance 
%figure('position',[672 695 1000 200]);
figure('position',[672 695 490 180]);
cc = get(groot,'defaultAxesColorOrder');
mrk = {'-o','-s','-v','-d','-^'};
xx = {-.2,-.1,0,.1,.2};
%subplot(2,1,1)

yyaxis left, hold on
decoder=1,
plot(X{decoder},SR{decoder},mrk{decoder},'markeredgecolor',cc(1,:),...
        'markersize',10,'markerfacecolor','w')


hold on 
for decoder=2:5,
    plot(X{decoder}(1),SR{decoder}(1),mrk{decoder},'markeredgecolor','w',...
        'markersize',8,'markerfacecolor',cc(1,:))
end
ylabel('% Success','FontSize',14)
%title('BCI Performance')
xlim([0.5,5.5])
ylim([0,1.1])
set(gca,'yTick',0:0.5:1,'yTickLabel',{'0' '50' '100'},'FontSize',16)
set(gca,'XTick',1:5,'XTickLabel',{'AM_1','PM_1','AM_2','PM_2','AM_3'},'FontSize',16)

% Stat check to see the data is within 2std
check=0;
for decoder=2:5
    if SR{decoder}(1)>mean(SR{1})-2*std(SR{1}) && SR{decoder}(1)<mean(SR{1})+2*std(SR{1});
        check=check+1;
    end
end

yyaxis right, hold on
decoder=1,
TtT = mean([TtST{decoder};TtRT{decoder}]);
plot(X{decoder},TtT,mrk{decoder},'markeredgecolor',cc(2,:),...
        'markersize',10,'markerfacecolor','w')

for decoder=2:5,
    TtT = mean([TtST{decoder};TtRT{decoder}]);
    plot(X{decoder}(1),mean(TtT),mrk{decoder},'markeredgecolor','w',...
        'markersize',8,'markerfacecolor',cc(2,:))
end
ylabel('Time to Target (sec)','FontSize',12)
%title('BCI Performance')
xlim([0.5,5.5])
ylim([0,15])
%set(gca,'XTick',0:5,'XTickLabel',{'AM_1','PM_1','AM_2','PM_2','AM_3'},'Clipping','off')
xlabel('Day','FontSize',16)

% Stat check to see the data is within 2std
check=0;
DecoderA=mean([TtST{1};TtRT{1}]);
for decoder=2:5
    TtT = mean([TtST{decoder};TtRT{decoder}]);
    if mean(TtT)>mean(DecoderA)-2*std(DecoderA) && mean(TtT)<mean(DecoderA)+2*std(DecoderA)
        check=check+1;
    end
end

HighQualityFigs('BCIPeformance')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[672 695 490 180]);
cc = get(groot,'defaultAxesColorOrder');
mrk = {'-o','-s','-v','-d','-^'};
xx = {-.2,-.1,0,.1,.2};
%subplot(2,1,2), hold on
%yyaxis
cla reset

decoder=1;
    plot(X{decoder},ITR{decoder},mrk{decoder},...
        'color','k','markeredgecolor','k','markersize',10,'markerfacecolor','w');

hold on;

for decoder=2:5,
    plot(X{decoder}(1),ITR{decoder}(1),mrk{decoder},...
        'color','k','markeredgecolor','w','markersize',8,'markerfacecolor','k');
end
%ylabel('Bits / sec','FontSize',16)
%xlabel('Day','FontSize',16)
%title('Fitt''s ITR')
xlim([0.5,5.5])
ylim([0.1,0.7])
set(gca,'yTick',[0.2,0.4,0.6],'yTickLabel',{'0.2', '0.4','0.6'},'FontSize',16,'box','off')
set(gca,'XTick',1:5,'XTickLabel',{'AM_1','PM_1','AM_2','PM_2','AM_3'},'FontSize',16,'box','off')
x = [0.5 5.5 5.5 0.5];
y = [0.3958-(2*0.0489) 0.3958-(2*0.0489) 0.3958+(2*0.0489) 0.3958+(2*0.0489)];
patch(x,y,'b','FaceAlpha',.3,'EdgeAlpha',.01)

HighQualityFigs('FittsRate')

% performing ttest between data points
[h,p,ci,stats]=ttest2(ITR{1},[ITR{2}(1),ITR{3}(1),ITR{4}(1),ITR{5}(1)])

%set(gca,'XTick',1:4,'XTickLabel',0:.5:4,'YTick',[.1,.55])
%legend(h,{'A','B','C','D','E'},'Position',[0.8746 0.1590 0.1153 0.2647])

% % save fig
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/PlugAndPlay/PerformancePerSessPerDecoderV2.pdf',...
%     '-pdf','-png','-r300','-painters')
% close(gcf)



%% does decoder move further away w/ sessions - timeseries angle
expts = [];

expts(end+1).yymmdd = '20190724';
expts(end).hhmmss = {'134612'};
expts(end).controller = {1};

expts(end+1).yymmdd = '20190724';
expts(end).hhmmss = {'140804'};
expts(end).controller = {2};

expts(end+1).yymmdd = '20190725';
expts(end).hhmmss = {'112929'};
expts(end).controller = {3};

expts(end+1).yymmdd = '20190725';
expts(end).hhmmss = {'140548'};
expts(end).controller = {4};

expts(end+1).yymmdd = '20190726';
expts(end).hhmmss = {'113418'};
expts(end).controller = {5};

% expts(end+1).yymmdd = '20190730';
% expts(end).hhmmss = {'112918'};
% expts(end).controller = {6};

% collect C matrix
C = cell(1,4);

% go through expts
ct = 0;
success = 1;
success_prev = 1;
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    for ii=1:length(expt.hhmmss),
        hhmmss = expt.hhmmss{ii};
        ct = ct + 1;
        
        % go through datafiles in fixed blocks
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        fprintf('\n%s-%s\n',yymmdd,hhmmss)
        files = dir(fullfile(datadir,'Data*.mat'));
        iii = length(files);
            
        % load data
        load(fullfile(datadir,files(iii).name));
        C{ct} = TrialData.KalmanFilter{1}.C;
    end
end

ref_idx = 1;
ang_x = [];
ang_y = [];
ang_c = [];
for i=1:5,
    ang_x(i) = subspace(C{ref_idx}(:,3),C{i}(:,3))*180/pi;
    ang_y(i) = subspace(C{ref_idx}(:,4),C{i}(:,4))*180/pi;
    ang_c(i) = subspace(C{ref_idx}(:,5),C{i}(:,5))*180/pi;
end
X = [1,2,3,4];
XX = [0,1,2,3,4];

% stats regressions
[h1,p1] = ttest(ang_x(2:end)');
mu1 = mean(ang_x(2:end));
[b1,bint1,~,~,S1] = regress(ang_x(2:end)',[ones(4,1),X']);
fprintf('\nVx:\n')
fprintf('t-test - eff_sz=%.03f,   p=%.03f\n',mu1,p1)
fprintf('regression - int=%.03f, slope=%.03f, CI: [%.03f, %.03f]\n',...
    b1(1),b1(2),bint1(2,1),bint1(2,2))

[h2,p2] = ttest(ang_y(2:end)');
mu2 = mean(ang_y(2:end));
[b2,bint2,~,~,S2] = regress(ang_y(2:end)',[ones(4,1),X']);
fprintf('\nVy:\n')
fprintf('regression - int=%.03f, slope=%.03f, CI: [%.03f, %.03f]\n',...
    b2(1),b2(2),bint2(2,1),bint2(2,2))

[h3,p3] = ttest(ang_c(2:end)');
mu3 = mean(ang_c(2:end));
[b3,bint3,~,~,S3] = regress(ang_c(2:end)',[ones(4,1),X']);
fprintf('\nConst:\n')
fprintf('regression - int=%.03f, slope=%.03f, CI: [%.03f, %.03f]\n',...
    b3(1),b3(2),bint3(2,1),bint3(2,2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
figure('position',[672 695 490 180]);
cc = get(groot,'defaultAxesColorOrder');
mrk = {'-o','-s','-v','-d','-^'};
xx = {-.2,-.1,0,.1,.2};

%subplot(2,1,1); hold on
% plot([1,4]-.1,b1(1)+b1(2)*[1,4],'-','color',cc(4,:))
decoder=1,
    plot(XX(decoder),ang_x(decoder),mrk{decoder},'markeredgecolor',cc(4,:),...
        'markersize',10,'markerfacecolor','w')
hold on
for decoder=2:5,
    plot(XX(decoder),ang_x(decoder),mrk{decoder},'markeredgecolor','w',...
        'markersize',10,'markerfacecolor',cc(4,:))
    
hold on
end
plot(XX,ang_x,'-','color',cc(4,:))
xlim([-.1,4.5])
ylim([-.5,6])

hold on;
% plot([1,4],b2(1)+b2(2)*[1,4],'-','color',cc(5,:))
decoder=1,
plot(XX(decoder),ang_y(decoder),mrk{decoder},'markeredgecolor',cc(5,:),...
        'markersize',10,'markerfacecolor','w')

hold on;
for decoder=2:5,
    plot(XX(decoder),ang_y(decoder),mrk{decoder},'markeredgecolor','w',...
        'markersize',10,'markerfacecolor',cc(5,:))
end
plot(XX,ang_y,'-','color',cc(5,:))

xlim([-.1,4.5])
ylim([-.5,6])
set(gca,'yTick',[0,6],'yTickLabel',{['0',char(176)],['6',char(176)]},'FontSize',16,'box','off')
set(gca,'XTick',0:4,'XTickLabel',{'A','B','C','D','E'},'FontSize',16,'box','off')
xlabel('Decoder','FontSize',16)
ylabel('\angle Decoder wrt A','FontSize',14)

HighQualityFigs('DecodersAngles')




% plot([1,4]+.1,b3(1)+b3(2)*[1,4],'-','color',cc(6,:))
% for decoder=1:5,
%     plot(XX(decoder),ang_c(decoder),mrk{decoder},'markeredgecolor','w',...
%         'markersize',6,'markerfacecolor',cc(6,:))
% end
% plot(XX,ang_c,'-','color',cc(6,:))

xlabel('Decoder')
ylabel('\angle KF wrt KF_1 (^{\circ})')
title('Decoder Stability Across Sessions')
set(gca,'XTick',0:4,'xticklabel',1:5,'XLim',[-0.5,4.5],'YLim',[0,6],'YTick',[0,6])

% save fig
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/PlugAndPlay/DecoderWeightsDivergence.pdf',...
%     '-pdf','-png','-r300','-painters')
% close(gcf)
