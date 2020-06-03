%% click trajectory analysis
clear, clc, close all
cc = get(groot,'defaultAxesColorOrder');

%% experiment info
expts = [];

expts(end+1).yymmdd = '20190819'; % best so far: being used in paper
expts(end).hhmmss = {'115346'};

% expts(end+1).yymmdd = '20190830'; % best so far: being used in paper
% expts(end).hhmmss = {'143231'};

%% saving
savedir = 'Figures/ClickerTrajAnalysis';

%% load clicker
% f = load('clicker_svm_mdl_812.mat');
f = load('clicker_svm_mdl_819.mat');
% f = load('clicker_svm_mdl_823.mat');
clicker = f.model;

%% plot example trajectory in click space

expt = expts(1);
yymmdd = expt.yymmdd;
hhmmss = expt.hhmmss{1};
        
% go through datafiles in fixed blocks
datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
    'GangulyServer','CursorControlGridTask',yymmdd,hhmmss,'BCI_Fixed');
fprintf('\n%s-%s\n',yymmdd,hhmmss)
files = dir(fullfile(datadir,'Data*.mat'));
        
% load data
trial = 1;
load(fullfile(datadir,files(trial).name));
time = TrialData.Time - TrialData.Time(1);

% run features through clicker
features = cat(2,TrialData.NeuralFeatures{:});
features = features(129:end,:);

thresh = TrialData.Params.DecisionBoundary;
traj = clicker.w * features - thresh;

mae = mean((traj(1:end-1)));

figure('units','normalized','position',[0.3516 0.7157 0.2333 0.1852])
hold on
xlim([time(1)-0.5,time(end)+0.5])
% ylim([-0.5,1.5])

stem(time(1:end-1), traj(1:end-1))
stem(time(end), traj(end))
plot([time(1),time(end-1)], repmat(mean(traj),1,2), '-', ...
    'color', [.6,.6,.6])
plot([time(1),time(end-1)], repmat(min(traj(1:end-1)),1,2), '-', ...
    'color', [.6,.6,.6])
plot([time(end),time(end)+0.5], repmat(min(traj),1,2), '-', ...
    'color', [.6,.6,.6])

title(sprintf('Neural Trajectory in Click Space: Trial %i', trial))
xlabel('Time (s)')
ylabel('Neural State')
text(0.5,0,'Click Threshold',...
    'verticalalignment','top')
text(time(end),mean(traj),'Mean',...
    'verticalalignment','middle')
text(time(end),min(traj(1:end-1)),'Min1',...
    'verticalalignment','middle')
text(time(end-1),min(traj),'Min2',...
    'verticalalignment','middle','horizontalalignment','right')

saveas(gcf, fullfile(savedir,'ExampleTrajectory'),'png')
export_fig(fullfile(savedir,'ExampleTrajectory'),...
    '-pdf','-r300','-painters')

%% repeat, looking at all trials in day, align to end of trial
expts = [];

expts(end+1).yymmdd = '20190819'; % best so far: being used in paper
expts(end).hhmmss = {'114110','114318','114616','115111','115346',...
    '142626','143335','143651','144113','144505','144816'};

% expts(end+1).yymmdd = '20190830'; % best so far: being used in paper
% expts(end).hhmmss = {'114031','115045','115346','115737','141547','141814',...
%     '142332','142709','143231','143449','143810','144326','144640','145048',...
%     '145507','145735','150124','150555'};

expt = expts(1);
yymmdd = expt.yymmdd;

win = 5;
trajs = [];
maes = [];
mns = [];
for ii=1:length(expt.hhmmss),
    hhmmss = expt.hhmmss{ii};
    
    datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
        'GangulyServer','CursorControlGridTask',yymmdd,hhmmss,'BCI_Fixed');
    datafiles = dir(fullfile(datadir,'Data*.mat'));
    T = length(datafiles);
    
    for trial=1:T,
        % load data
        load(fullfile(datadir,datafiles(trial).name));
        time = TrialData.Time - TrialData.Time(1);
        
        features = cat(2,TrialData.NeuralFeatures{:});
        features = features(129:end,:);
        
        traj = clicker.w * features;
        thresh = TrialData.Params.DecisionBoundary;
        
        mae = mean((traj(1:end-1)-thresh));
        mn = min((traj(1:end-1)-thresh));
        
        % track
        if length(time) >= win*TrialData.Params.ScreenRefreshRate,
            tidx = time >= (time(end) - win);
            trajs(trial,:) = traj(tidx);
        else, % zero pad
            N = 5*TrialData.Params.ScreenRefreshRate - length(time);
            trajs(trial,:) = padarray(traj,[0,N],nan,'pre');
        end
        maes(trial,1) = mae;
        mns(trial,1) = mn;
    end
end

figure('units','normalized')
hold on
xlim([-5,0.2])
ylim([-1,1])

shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs,1),nanstd(trajs),{'color',cc(1,:)});
hline(thresh, 'r--')

title(sprintf('Average Neural Trajectory in Click Space (1 day)'))
xlabel('Time (s)')
ylabel('Neural State')
text(-win+0.5,thresh+0.1,'Click Threshold',...
    'verticalalignment','top')
text(-win+0.5,thresh-0.1,sprintf('Avg Dist from Thresh: %0.2f +/- %0.2fs.d.', mean(maes), std(maes)),...
    'verticalalignment','top')
text(-win+0.5,thresh-0.2,sprintf('*shaded region is s.d.'),...
    'verticalalignment','top')

% saveas(gcf, fullfile(savedir,'AvgTrajectoriesDay'),'png')

% continue trials plot
figure('position',[441 293 179 577]);

subplot(3,1,1)
plot(maes,'.')
hline(-0.5,'r--')
xlim([0,length(mns)+1])
ylim([-1,1])
title('avg dist per trial')

subplot(3,1,2)
plot(mns,'.')
hline(-0.5,'r--')
xlim([0,length(mns)+1])
ylim([-1,1])
title('min dist per trial')

subplot(3,1,3)
shadedErrorBar(-4.8:0.2:0,nanmean(trajs),nanstd(trajs),{'color',cc(1,:)})
hline(-0.5,'r--')
xlim([-5,0.2])
ylim([-1,1])
title('avg trajectory')
% saveas(gcf, fullfile(savedir,'ClickerDayAnalysis'),'png')


%% go through each day (before clicker was introduced) and collect clicking traj & dist

% experiment info
expts = [];

expts(end+1).yymmdd = '20190501';
expts(end).hhmmss = {'105345','112049','133745'};
expts(end).perf_hhmmss = {'105345','112049','135512'};

expts(end+1).yymmdd = '20190507';
expts(end).hhmmss = {'105331','112817'};
expts(end).perf_hhmmss = {'111526','113944'};

expts(end+1).yymmdd = '20190510';
expts(end).hhmmss = {'112633','140749'};
expts(end).perf_hhmmss = {'141900'};

expts(end+1).yymmdd = '20190514';
expts(end).hhmmss = {'112836','133050'};
expts(end).perf_hhmmss = {'135927','140808','141731','142704'};

expts(end+1).yymmdd = '20190515';
expts(end).hhmmss = {'110735'};
expts(end).perf_hhmmss = {'112302','113447'};

expts(end+1).yymmdd = '20190521';
expts(end).hhmmss = {'135008','141438'};
expts(end).perf_hhmmss = {'135943','141438'};

expts(end+1).yymmdd = '20190524';
expts(end).hhmmss = {'110219','133313'};
expts(end).perf_hhmmss = {'111731','112353','134653','135957'};

expts(end+1).yymmdd = '20190529';
expts(end).hhmmss = {'105609','113247'};
expts(end).perf_hhmmss = {'111528','114050'};

expts(end+1).yymmdd = '20190531';
expts(end).hhmmss = {'102944','111946','132244','135046'};
expts(end).perf_hhmmss = {'105048','110703','112444','133517','140204','141319'};

expts(end+1).yymmdd = '20190604';
expts(end).hhmmss = {'112454','140706'};
expts(end).perf_hhmmss = {'112454','114636','143109'};

expts(end+1).yymmdd = '20190607';
expts(end).hhmmss = {'104622','133540','140114'};
expts(end).perf_hhmmss = {'105615','135050','140828'};

expts(end+1).yymmdd = '20190610';
expts(end).hhmmss = {'104104'};
expts(end).perf_hhmmss = {'105809','110746'};

expts(end+1).yymmdd = '20190618';
expts(end).hhmmss = {'132907','135020','141917'};
expts(end).perf_hhmmss = {'133933','135944','143103'};

expts(end+1).yymmdd = '20190621';
expts(end).hhmmss = {'110109','112545','133227','135917'};
expts(end).perf_hhmmss = {'110942','113715','134143','141129'};

expts(end+1).yymmdd = '20190626';
expts(end).hhmmss = {'104807','105539','110553','111358','112324'};
expts(end).perf_hhmmss = {'105539','110229','110553','111358','112324','113744','133301','134240','135658'};

days = {expts.yymmdd};

% load data

% go through expts
trajs = {};
maes = {};
mns = {};
trial = 1;
trials = [];
day_break = [];
session_break = [];
session_init = [];
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    day_break(end+1) = trial; 
    
    for ii=1:length(expt.hhmmss),
        hhmmss = expt.hhmmss{ii};
        session_break(end+1) = trial; 
        
        datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_CLDA');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        T = length(datafiles);

        for iii=1:T,
            % load data, grab neural features
            disp(datafiles(iii).name)
            load(fullfile(datadir,datafiles(iii).name))
            if iii==1,
                session_init(end+1) = TrialData.Params.InitializationMode; 
            end
            
            % compute
            time = TrialData.Time - TrialData.Time(1);
            features = cat(2,TrialData.NeuralFeatures{:});
            features = features(129:end,:);

            traj = clicker.w * features;
            thresh = -0.5;

            mae = mean((traj(1:end-1)-thresh));
            mn = min((traj(1:end-1)-thresh));

            % align to click / end of trial
            if length(time) >= win*TrialData.Params.ScreenRefreshRate,
                tidx = time >= (time(end) - win);
                traj = traj(tidx);
            else, % zero pad
                N = win*TrialData.Params.ScreenRefreshRate - length(time);
                traj = padarray(traj,[0,N],nan,'pre');
            end
            
            % sample at 5hz
            if TrialData.Params.ScreenRefreshRate ~= 5,
                traj = downsample(traj,TrialData.Params.ScreenRefreshRate/5);
            end

            if length(traj)~=25,
            end
            
            % track
            trajs{trial} = traj;
            maes{trial} = mae;
            mns{trial} = mn;
            
            trials(trial) = trial;
            trial = trial + 1;
        end
        
    end
end
day_break(end+1) = trial;
session_break(end+1) = trial; 

% plot
trajs_mat = cat(1,trajs{:});
maes_mat = cat(1,maes{:});
mns_mat = cat(1,mns{:});

figure('position',[445 326 1088 547]);

% plot maes across days
subplot(3,1,1)
scatter(trials,maes_mat,'.')
hline(-0.5,'r--')
xlim([trials(1)-1,trials(end)+1])
ylim([-1,1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
title('avg dist per trial')

% plot mns across days
subplot(3,1,2)
scatter(trials,mns_mat,'.')
hline(-0.5,'r--')
xlim([trials(1)-1,trials(end)+1])
ylim([-1,1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
title('min dist per trial')

% plot ex trajs across days
subplot(3,3,7)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(1)) & (trials < day_break(2));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:),1),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 1')

% plot ex trajs across days
subplot(3,3,8)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(7)) & (trials < day_break(8));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:)),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 7')

% plot ex trajs across days
subplot(3,3,9)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(14)) & (trials < day_break(15));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:),1),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 15')

% saveas(gcf, fullfile(savedir,'ClickerAnalysisBeforeClicking_CLDA'),'png')

%% repeat for fixed blocks

% experiment info
expts = [];

expts(end+1).yymmdd = '20190501';
expts(end).hhmmss = {'105345','112049','133745'};
expts(end).perf_hhmmss = {'105345','112049','135512'};

expts(end+1).yymmdd = '20190507';
expts(end).hhmmss = {'105331','112817'};
expts(end).perf_hhmmss = {'111526','113944'};

expts(end+1).yymmdd = '20190510';
expts(end).hhmmss = {'112633','140749'};
expts(end).perf_hhmmss = {'141900'};

expts(end+1).yymmdd = '20190514';
expts(end).hhmmss = {'112836','133050'};
expts(end).perf_hhmmss = {'135927','140808','141731','142704'};

expts(end+1).yymmdd = '20190515';
expts(end).hhmmss = {'110735'};
expts(end).perf_hhmmss = {'112302','113447'};

expts(end+1).yymmdd = '20190521';
expts(end).hhmmss = {'135008','141438'};
expts(end).perf_hhmmss = {'135943','141438'};

expts(end+1).yymmdd = '20190524';
expts(end).hhmmss = {'110219','133313'};
expts(end).perf_hhmmss = {'111731','112353','134653','135957'};

expts(end+1).yymmdd = '20190529';
expts(end).hhmmss = {'105609','113247'};
expts(end).perf_hhmmss = {'111528','114050'};

expts(end+1).yymmdd = '20190531';
expts(end).hhmmss = {'102944','111946','132244','135046'};
expts(end).perf_hhmmss = {'105048','110703','112444','133517','140204','141319'};

expts(end+1).yymmdd = '20190604';
expts(end).hhmmss = {'112454','140706'};
expts(end).perf_hhmmss = {'112454','114636','143109'};

expts(end+1).yymmdd = '20190607';
expts(end).hhmmss = {'104622','133540','140114'};
expts(end).perf_hhmmss = {'105615','135050','140828'};

expts(end+1).yymmdd = '20190610';
expts(end).hhmmss = {'104104'};
expts(end).perf_hhmmss = {'105809','110746'};

expts(end+1).yymmdd = '20190618';
expts(end).hhmmss = {'132907','135020','141917'};
expts(end).perf_hhmmss = {'133933','135944','143103'};

expts(end+1).yymmdd = '20190621';
expts(end).hhmmss = {'110109','112545','133227','135917'};
expts(end).perf_hhmmss = {'110942','113715','134143','141129'};

expts(end+1).yymmdd = '20190626';
expts(end).hhmmss = {'104807','105539','110553','111358','112324'};
expts(end).perf_hhmmss = {'105539','110229','110553','111358','112324','113744','133301','134240','135658'};

days = {expts.yymmdd};

% load data

% go through expts
trajs = {};
maes = {};
mns = {};
trial = 1;
trials = [];
day_break = [];
session_break = [];
session_init = [];
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    day_break(end+1) = trial; 
    
    for ii=1:length(expt.perf_hhmmss),
        hhmmss = expt.perf_hhmmss{ii};
        session_break(end+1) = trial; 
        
        datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        T = length(datafiles);

        for iii=1:T,
            % load data, grab neural features
            disp(datafiles(iii).name)
            load(fullfile(datadir,datafiles(iii).name))
            if iii==1,
                session_init(end+1) = TrialData.Params.InitializationMode; 
            end
            
            % compute
            time = TrialData.Time - TrialData.Time(1);
            features = cat(2,TrialData.NeuralFeatures{:});
            features = features(129:end,:);

            traj = clicker.w * features;
            thresh = -0.5;

            mae = mean((traj(1:end-1)-thresh));
            mn = min((traj(1:end-1)-thresh));

            % align to click / end of trial
            if length(time) >= win*TrialData.Params.ScreenRefreshRate,
                tidx = time >= (time(end) - win);
                traj = traj(tidx);
            else, % zero pad
                N = win*TrialData.Params.ScreenRefreshRate - length(time);
                traj = padarray(traj,[0,N],nan,'pre');
            end
            
            % sample at 5hz
            if TrialData.Params.ScreenRefreshRate ~= 5,
                traj = downsample(traj,TrialData.Params.ScreenRefreshRate/5);
            end

            if length(traj)~=25,
            end
            
            % track
            trajs{trial} = traj;
            maes{trial} = mae;
            mns{trial} = mn;
            
            trials(trial) = trial;
            trial = trial + 1;
        end
        
    end
end
day_break(end+1) = trial;
session_break(end+1) = trial; 

% plot
trajs_mat = cat(1,trajs{:});
maes_mat = cat(1,maes{:});
mns_mat = cat(1,mns{:});

figure('position',[445 326 1088 547]);

% plot maes across days
subplot(3,1,1)
scatter(trials,maes_mat,'.')
hline(-0.5,'r--')
xlim([trials(1)-1,trials(end)+1])
ylim([-1,1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
title('avg dist per trial')

% plot mns across days
subplot(3,1,2)
scatter(trials,mns_mat,'.')
hline(-0.5,'r--')
xlim([trials(1)-1,trials(end)+1])
ylim([-1,1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
title('min dist per trial')

% plot ex trajs across days
subplot(3,3,7)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(1)) & (trials < day_break(2));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:),1),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 1')

% plot ex trajs across days
subplot(3,3,8)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(7)) & (trials < day_break(8));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:)),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 7')

% plot ex trajs across days
subplot(3,3,9)
hold on
xlim([-5,0.2])
ylim([-1,1])
idx = (trials >= day_break(14)) & (trials < day_break(15));
shadedErrorBar(-(win-0.2):0.2:0,nanmean(trajs_mat(idx,:),1),nanstd(trajs_mat(idx,:)),{'color',cc(1,:)});
hline(thresh, 'r--')
title('day 15')

saveas(gcf, fullfile(savedir,'ClickerAnalysisBeforeClicking_Fixed'),'png')

%% pre-clicking final day vs clicking
% experiment info
expts = [];

expts(end+1).yymmdd = '20190626';
expts(end).hhmmss = {'104807','105539','110553','111358','112324'};
expts(end).perf_hhmmss = {'105539','110229','110553','111358','112324','113744','133301','134240','135658'};

% go through expts
trajs = {};
maes = {};
mns = {};
trial = 1;
trials = [];
day_break = [];
session_break = [];
session_init = [];
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    day_break(end+1) = trial; 
    
    for ii=1:length(expt.perf_hhmmss),
        hhmmss = expt.perf_hhmmss{ii};
        session_break(end+1) = trial; 
        
        datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        T = length(datafiles);

        for iii=1:T,
            % load data, grab neural features
            disp(datafiles(iii).name)
            load(fullfile(datadir,datafiles(iii).name))
            if iii==1,
                session_init(end+1) = TrialData.Params.InitializationMode; 
            end
            
            % compute
            time = TrialData.Time - TrialData.Time(1);
            features = cat(2,TrialData.NeuralFeatures{:});
            features = features(129:end,:);

            traj = clicker.w * features;
            thresh = -0.5;

            mae = mean((traj(1:end-1)-thresh));
            mn = min((traj(1:end-1)-thresh));

            % align to click / end of trial
            if length(time) >= win*TrialData.Params.ScreenRefreshRate,
                tidx = time >= (time(end) - win);
                traj = traj(tidx);
            else, % zero pad
                N = win*TrialData.Params.ScreenRefreshRate - length(time);
                traj = padarray(traj,[0,N],nan,'pre');
            end
            
            % sample at 5hz
            if TrialData.Params.ScreenRefreshRate ~= 5,
                traj = downsample(traj,TrialData.Params.ScreenRefreshRate/5);
            end

            if length(traj)~=25,
            end
            
            % track
            trajs{trial} = traj;
            maes{trial} = mae;
            mns{trial} = mn;
            
            trials(trial) = trial;
            trial = trial + 1;
        end
        
    end
end
day_break(end+1) = trial;
session_break(end+1) = trial; 

pre_trajs = cat(1,trajs{:});
pre_maes = cat(1,maes{:});
pre_mns = cat(1,mns{:});

expts = [];

expts(end+1).yymmdd = '20190819'; % best so far: being used in paper
expts(end).hhmmss = {'114110','114318','114616','115111','115346',...
    '142626','143335','143651','144113','144505','144816'};

% expts(end+1).yymmdd = '20190830'; % best so far: being used in paper
% expts(end).hhmmss = {'114031','115045','115346','115737','141547','141814',...
%     '142332','142709','143231','143449','143810','144326','144640','145048',...
%     '145507','145735','150124','150555'};

expt = expts(1);
yymmdd = expt.yymmdd;

win = 5;
trajs = [];
post_maes = [];
post_mns = [];
for ii=1:length(expt.hhmmss),
    hhmmss = expt.hhmmss{ii};
    
    datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
        'GangulyServer','CursorControlGridTask',yymmdd,hhmmss,'BCI_Fixed');
    datafiles = dir(fullfile(datadir,'Data*.mat'));
    T = length(datafiles);
    
    for trial=1:T,
        % load data
        load(fullfile(datadir,datafiles(trial).name));
        time = TrialData.Time - TrialData.Time(1);
        
        features = cat(2,TrialData.NeuralFeatures{:});
        features = features(129:end,:);
        
        traj = clicker.w * features;
        thresh = TrialData.Params.DecisionBoundary;
        
        mae = mean((traj(1:end-1)-thresh));
        mn = min((traj(1:end-1)-thresh));
        
        % track
        if length(time) >= win*TrialData.Params.ScreenRefreshRate,
            tidx = time >= (time(end) - win);
            post_trajs(trial,:) = traj(tidx);
        else, % zero pad
            N = 5*TrialData.Params.ScreenRefreshRate - length(time);
            post_trajs(trial,:) = padarray(traj,[0,N],nan,'pre');
        end
        post_maes(trial,1) = mae;
        post_mns(trial,1) = mn;
    end
end

% plots
figure()

subplot(2,2,1), hold on
histogram(pre_maes,10,'normalization','pdf')
histogram(post_maes,10,'normalization','pdf')
ylabel('pdf')
legend('Pre-Clicking','Clicking')
title('avg dist from thresh')

subplot(2,2,2), hold on
histogram(pre_mns,10,'normalization','pdf')
histogram(post_mns,10,'normalization','pdf')
ylabel('pdf')
legend('Pre-Clicking','Clicking')
title('min point in click space')

subplot(2,2,3), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(pre_trajs),nanstd(pre_trajs),{'color',cc(1,:)})
hline(-0.5,'r--')
ylim([-1,1])
xlim([-5,0.2])
ylabel('Clicker Space')
title('Pre-Clicking Trajectories')

subplot(2,2,4), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(post_trajs),nanstd(post_trajs),{'color',cc(1,:)})
hline(-0.5,'r--')
ylim([-1,1])
xlim([-5,0.2])
ylabel('Clicker Space')
title('Clicking Trajectories')

saveas(gcf, fullfile(savedir,'PreVsPostClickerComparison'),'png')

%% pre-clicking 4 days vs clicking
% experiment info
expts = [];

expts(end+1).yymmdd = '20190618';
expts(end).hhmmss = {'132907','135020','141917'};
expts(end).perf_hhmmss = {'133933','135944','143103'};

expts(end+1).yymmdd = '20190621';
expts(end).hhmmss = {'110109','112545','133227','135917'};
expts(end).perf_hhmmss = {'110942','113715','134143','141129'};

expts(end+1).yymmdd = '20190626';
expts(end).hhmmss = {'104807','105539','110553','111358','112324'};
expts(end).perf_hhmmss = {'105539','110229','110553','111358','112324','113744','133301','134240','135658'};

win = 5;

% go through expts
trajs = {};
maes = {};
mns = {};
trial = 1;
trials = [];
day_break = [];
session_break = [];
session_init = [];
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    day_break(end+1) = trial; 
    
    for ii=1:length(expt.perf_hhmmss),
        hhmmss = expt.perf_hhmmss{ii};
        session_break(end+1) = trial; 
        
        datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,hhmmss,'BCI_Fixed');
        datafiles = dir(fullfile(datadir,'Data*.mat'));
        T = length(datafiles);

        for iii=1:T,
            % load data, grab neural features
            disp(datafiles(iii).name)
            load(fullfile(datadir,datafiles(iii).name))
            if iii==1,
                session_init(end+1) = TrialData.Params.InitializationMode; 
            end
            
            % compute
            time = TrialData.Time - TrialData.Time(1);
            features = cat(2,TrialData.NeuralFeatures{:});
            features = features(129:end,:);

            thresh = -0.5;
            traj = (clicker.w * features) - thresh;

            mae = mean(traj);
            mn = min(traj);

            % align to click / end of trial
            if length(time) >= win*TrialData.Params.ScreenRefreshRate,
                tidx = time >= (time(end) - win);
                traj = traj(tidx);
            else, % zero pad
                N = win*TrialData.Params.ScreenRefreshRate - length(time);
                traj = padarray(traj,[0,N],nan,'pre');
            end
            
            % sample at 5hz
            if TrialData.Params.ScreenRefreshRate ~= 5,
                traj = downsample(traj,TrialData.Params.ScreenRefreshRate/5);
            end

            if length(traj)~=25,
            end
            
            % track
            trajs{trial} = traj;
            maes{trial} = mae;
            mns{trial} = mn;
            
            trials(trial) = trial;
            trial = trial + 1;
        end
        
    end
end
day_break(end+1) = trial;
session_break(end+1) = trial; 

pre_trajs = cat(1,trajs{:});
pre_maes = cat(1,maes{:});
pre_mns = cat(1,mns{:});

expts = [];

expts(end+1).yymmdd = '20190819'; % best so far: being used in paper
expts(end).hhmmss = {'114110','114318','114616','115111','115346',...
    '142626','143335','143651','144113','144505','144816'};

% expts(end+1).yymmdd = '20190830'; % best so far: being used in paper
% expts(end).hhmmss = {'114031','115045','115346','115737','141547','141814',...
%     '142332','142709','143231','143449','143810','144326','144640','145048',...
%     '145507','145735','150124','150555'};

expt = expts(1);
yymmdd = expt.yymmdd;

win = 5;
trajs = [];
post_maes = [];
post_mns = [];
post_mns2 = [];
for ii=1:length(expt.hhmmss),
    hhmmss = expt.hhmmss{ii};
    
    datadir = fullfile('/media/dsilver/data/Bravo1',yymmdd,...
        'GangulyServer','CursorControlGridTask',yymmdd,hhmmss,'BCI_Fixed');
    datafiles = dir(fullfile(datadir,'Data*.mat'));
    T = length(datafiles);
    
    for trial=1:T,
        % load data
        load(fullfile(datadir,datafiles(trial).name));
        time = TrialData.Time - TrialData.Time(1);
        
        features = cat(2,TrialData.NeuralFeatures{:});
        features = features(129:end,:);
        
        thresh = TrialData.Params.DecisionBoundary;
        traj = (clicker.w * features) - thresh;

        mae = mean(traj(1:end-1));
        mn = min(traj(1:end-1));
        
        % track
        if length(time) >= win*TrialData.Params.ScreenRefreshRate,
            tidx = time >= (time(end) - win);
            post_trajs(trial,:) = traj(tidx);
        else, % zero pad
            N = 5*TrialData.Params.ScreenRefreshRate - length(time);
            post_trajs(trial,:) = padarray(traj,[0,N],nan,'pre');
        end
        post_maes(trial,1) = mae;
        post_mns(trial,1) = mn;
        post_mns2(trial,1) = min(traj);
    end
end

% plots
figure('position',[681 600 731 350])

subplot(2,1,1), hold on
plot(trials,pre_mns,'.','color',cc(1,:))
h(1)=plot(trials(end)+10+(1:length(post_mns)),post_mns,'.','color',cc(1,:));
h(2)=plot(trials(end)+10+(1:length(post_mns)),post_mns2,'.','color',cc(2,:));
vline([day_break(2:end),day_break(end)+10],'-k','linewidth',1.5)
ylabel('Min')
l=legend(h,{'Min1','Min2'},'location','southeast','position',[0.7136 0.6010 0.1094 0.0986]);
ylim([-1,0.5])
xlim([0,trials(end)+11+length(post_mns)])

idx = trials < day_break(2);
subplot(2,4,5), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(pre_trajs(idx,:)),...
    nanstd(pre_trajs(idx,:))./sqrt(sum(idx)),{'color',cc(1,:)})
hline(0,'r-')
ylim([-0.25,1.25])
xlim([-5,0.2])
ylabel('Clicker Space')
title('Day +28')

idx = trials < day_break(3) & trials >= day_break(2);
subplot(2,4,6), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(pre_trajs(idx,:)),...
    nanstd(pre_trajs(idx,:))./sqrt(sum(idx)),{'color',cc(1,:)})
hline(0,'r-')
ylim([-0.25,1.25])
xlim([-5,0.2])
title('Day +31')
set(gca,'YTick',[])

idx = trials < day_break(4) & trials >= day_break(3);
subplot(2,4,7), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(pre_trajs(idx,:)),...
    nanstd(pre_trajs(idx,:))./sqrt(sum(idx)),{'color',cc(1,:)})
hline(0,'r-')
ylim([-0.25,1.25])
xlim([-5,0.2])
title('Day +36')
set(gca,'YTick',[])

subplot(2,4,8), hold on
shadedErrorBar(-4.8:0.2:0,nanmean(post_trajs),...
    nanstd(post_trajs)./sqrt(size(post_trajs,1)),{'color',cc(1,:)})
hline(0,'r-')
ylim([-0.25,1.25])
xlim([-5,0.2])
title('Clicker Day')
set(gca,'YTick',[])


saveas(gcf, fullfile(savedir,'PreVsPostClickerComparison'),'png')
export_fig(fullfile(savedir,'PreVsPostClickerComparison'),...
    '-pdf','-r300','-painters')
