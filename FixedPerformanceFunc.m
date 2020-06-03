
function [FittsRate, SuccessRate, TravelTime,Clicker]=FixedPerformanceFunc(yymmdd,hhmmss,FixedType)


%% Calculate the performance in two tasks

if strcmp(FixedType,'Center-Out')
    
    target_ang = 0:45:360-45;
    
    % go through expts
    success = 1;
    success_prev = 1;
    
    % per session
    SR_sess     = cell(1,8); % Success Rate
    TT_sess     = cell(1,8); % Total Targets
    TtST_sess   = cell(1,8); % Time to Start Target
    TtRT_sess   = cell(1,8); % Time to Reach Target
    
    % go through datafiles in fixed blocks  
    fprintf('Fixed_Center-Out:')
    fprintf('\n%s-%s\n',yymmdd,hhmmss)
    
    if str2num(yymmdd)< 20191127 % the structure of folders changed
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer',FixedType,yymmdd,hhmmss,'BCI_Fixed');
    else
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer',yymmdd,'CenterOut',hhmmss,'BCI_Fixed'); 
    end
    
    files = dir(fullfile(datadir,'Data*.mat'));
    for iii=1:length(files),
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
    
    
    % Output
    %fprintf('Trials: %i\n',length(files))
    %fprintf('Target Radius: %i\n',TrialData.Params.ReachTargetRadius)
    %fprintf('Center-Reset: %i\n',TrialData.Params.CenterReset)
    %fprintf('A: %.03f\n',TrialData.Params.KF.A(3,3))
    %fprintf('W: %i\n',TrialData.Params.KF.W(3,3))
    %fprintf('G: %i\n',TrialData.Params.Gain)
    
    %fprintf('Success Rate: %i%%\n',round(100*AvgSR_sess))
    %fprintf('Time to Target: %.1f sec\n',nanmean([AvgTtST,AvgTtRT]))
    %fprintf('Fitt''s IT Rate: %.2f bits/sec\n',ITR_sess)
    
    FittsRate=ITR_sess;
    SuccessRate=round(100*AvgSR_sess);
    TravelTime=nanmean([AvgTtST,AvgTtRT]);
    Clicker='off';
    
    
elseif strcmp(FixedType,'CursorControlGridTask') || strcmp(FixedType,'GridTask')  
    
        % go through datafiles in fixed blocks  
    fprintf('Fixed_GridTask:')
    fprintf('\n%s-%s\n',yymmdd,hhmmss)
    
     datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
         'GangulyServer',FixedType,yymmdd,hhmmss,'BCI_Fixed');
    
    files = dir(fullfile(datadir,'Data*.mat'));
    load(fullfile(datadir,files(1).name));
    
    if isfield(TrialData,'ClickerState')
        Clicker='on';
        
        acc=0;fail=0;t=[];fail_index=[];
        for iii=1:length(files)
            load(fullfile(datadir,files(iii).name));
            
            if TrialData.TargetID ==  TrialData.SelectedTargetID
                acc=acc+1;
            elseif TrialData.TargetID ~=  TrialData.SelectedTargetID ...
                    && ~TrialData.SelectedTargetID == 0
                fail = fail+1;
                fail_index = [fail_index iii];
            end
            
            %t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
            t=[t (length(TrialData.Time)-1)*(1/TrialData.Params.UpdateRate)];
        end
        bitrate = (log2(prod(TrialData.Params.GridLayout))*max(0,acc-fail))/sum(t);
        FittsRate=bitrate;
        SuccessRate=acc/length(files)*100;
        TravelTime=nanmean(t);
        
    elseif ~isfield(TrialData,'ClickerState')
        Clicker='off';
        
        acc=0;t=[];
        for iii=1:length(files)
            load(fullfile(datadir,files(iii).name));
            
            if TrialData.ErrorID==0
                acc=acc+1;
            end
            
            t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
        end
        bitrate = (log2(prod(TrialData.Params.GridLayout))*max(0,acc))/sum(t);
        FittsRate=bitrate;
        SuccessRate=acc/length(files)*100;
        TravelTime=nanmean(t);
        
    end
    
    
elseif strcmp(FixedType,'RadialTask')
    
    % go through datafiles in fixed blocks
    fprintf('Fixed_RadialTask:')
    fprintf('\n%s-%s\n',yymmdd,hhmmss)
    
    datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
        'GangulyServer',yymmdd,FixedType,hhmmss,'BCI_Fixed');
    
    files = dir(fullfile(datadir,'Data*.mat'));
    load(fullfile(datadir,files(1).name));
    
    if isfield(TrialData,'ClickerState')
        Clicker=1;%'on';
        
        acc=0;fail=0;t=[];fail_index=[];correct_index=[];
        for i=1:length(files)
            load(fullfile(datadir,files(i).name));
            if (length(TrialData.Time)-1)*1/TrialData.Params.UpdateRate < 20
                if TrialData.TargetID ==  TrialData.SelectedTargetID
                    acc=acc+1;
                    correct_index = [correct_index i];
                elseif TrialData.TargetID ~=  TrialData.SelectedTargetID ...
                        && ~TrialData.SelectedTargetID == 0
                    fail = fail+1;
                    fail_index = [fail_index i];
                end
                %t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
                t=[t (length(TrialData.Time)-1)*1/TrialData.Params.UpdateRate];
            end
        end
        bitrate=(3*max(0,acc-fail))/sum(t);
        FittsRate=bitrate;
        SuccessRate=acc/length(files)*100;
        TravelTime=nanmean(t);
        
    elseif ~isfield(TrialData,'ClickerState')
        Clicker='off';
        
        acc=0;t=[];
        for i=1:length(files)
            load(fullfile(datadir,files(i).name));
            
            if TrialData.ErrorID==0
                acc=acc+1;
            end
            
            t=[t TrialData.TrialEndTime-TrialData.TrialStartTime];
        end
        bitrate =(3*max(0,acc-fail))/sum(t);
        FittsRate=bitrate;
        SuccessRate=acc/length(files)*100;
        TravelTime=nanmean(t);
        
    end
    
end


end


