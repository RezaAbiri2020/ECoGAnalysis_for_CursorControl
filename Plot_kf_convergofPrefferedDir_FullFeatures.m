
%% check if preffered direction of kalman weights are converging within and across days, compare
% with performance
% analysis of angle btw velocity vectors 
clear, clc, close all
cc = get(groot,'defaultAxesColorOrder');

%% experiment info
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

%% load performance metrics
fid = fopen('performance_v2.csv');
perf = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f',...
    'Delimiter',',','Headerlines',1);
fclose(fid);

perf_idx = [];
ct = 0;
for i=1:length(expts),
    expt = expts(i);
    yymmdd = expt.yymmdd;
    for ii=1:length(expt.perf_hhmmss),
        perf_hhmmss = str2num(expt.perf_hhmmss{ii});
        
        % find previous rml session
        hhmmss = cellfun(@str2num,expt.hhmmss,'un',0);
        hhmmss = [hhmmss{:}];
        idx = find(hhmmss<=perf_hhmmss,1,'last');
        
        % store
        if isempty(idx), % use prev
            perf_idx(end+1) = perf_idx(end);
        else,
            perf_idx(end+1) = ct + idx;
        end
        
    end
    % count of rml sessions
    ct = ct + length(hhmmss);
end

idx = setdiff(1:length(perf{1})-3,6); % ignore masks
perf_idx = 1 + perf_idx(idx);

% % store measures in cell array
% Perf{ct,1} = sprintf('%s-%s',yymmdd,hhmmss);
% Perf{ct,2} = length(files);
% Perf{ct,3} = TrialData.Params.KF.A(3,3);
% Perf{ct,4} = TrialData.Params.KF.W(3,3);
% Perf{ct,5} = TrialData.Params.Gain;
% Perf{ct,6} = TrialData.Params.TargetSize;
% Perf{ct,7} = TrialData.Params.ReachTargetRadius;
% Perf{ct,8} = TrialData.Params.CenterReset;
% Perf{ct,9} = TrialData.Params.CLDA.Lambda;
% Perf{ct,10} = TrialData.Params.CLDA.FinalLambda;
% Perf{ct,11} = expt.masks{ii};
% Perf{ct,12} = round(100*AvgSR_sess);
% Perf{ct,13} = round(100*AvgRR_sess);
% Perf{ct,14} = AvgTtST;
% Perf{ct,15} = AvgTtRT;
% Perf{ct,16} = nanmean([AvgTtST,AvgTtRT]);
% Perf{ct,17} = ITR_RT_sess;
% Perf{ct,18} = ITR_sess;
% Perf{ct,19} = AvgJerk_sess;

%% load data

% go through expts
Ccell = {};
alpha = {};
lambda = {};
TargetID = {};
TargetIDstr = {};
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
        
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data',yymmdd,...
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
            
            Ccell{trial} = TrialData.KalmanFilter{1}.C;
            alpha{trial} = TrialData.CursorAssist(1);
            %lambda{trial} = round(TrialData.KalmanFilter{1}.Lambda,2);
            TargetID{trial} = TrialData.TargetID;
            TargetIDstr{trial} = int2str(TrialData.TargetID);
            
            trials(trial) = trial;
            trial = trial + 1;
        end
        
    end
end
day_break(end+1) = trial;
session_break(end+1) = trial; 

%% plot behavioral performance
idx = setdiff(1:length(perf{1})-3,6); % ignore masks
sz = 30;

figure('position',[445 326 1088 547]);

subplot(4,1,1)
yyaxis left
hold on
ylim([0,100])
xlim([0,trials(end)+1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'-','color',[.7,.7,.7]);
end
scatter(session_break(perf_idx),perf{12}(idx),'o','markeredgecolor','w',...
    'sizedata',sz,'markerfacecolor',cc(1,:)); hold on
plot(session_break(perf_idx),smooth(session_break(perf_idx),perf{12}(idx),20),'--','color',cc(1,:))
ylabel('% Success')

yyaxis right
hold on
scatter(session_break(perf_idx),perf{15}(idx),'o','markeredgecolor','w',...
    'sizedata',sz,'markerfacecolor',cc(2,:)); hold on
plot(session_break(perf_idx),smooth(session_break(perf_idx),perf{15}(idx),20),'--','color',cc(2,:))
ylim([0,15])
xlim([0,trials(end)+1])
ylabel('Time to Target (sec)')
set(gca,'XTick',[])
title('BCI Performance')
set(gca,'Clipping','off')

% STATS

% t-test- FIRST 5 vs LAST 5 DAYS
ITR = perf{17}(idx);
ITR_daily = ITR(1:12);
ITR_cont = ITR(end-20:end);
[h1,p1] = ttest2(ITR_daily,ITR_cont);
mu1 = mean(ITR_cont) - mean(ITR_daily);

% regression stats - FIRST 5 vs MIDDLE 5 DAYS (not enough data pts to be
% significant. I'm ignoring even though slopes are prob diff)
ITR = perf{17}(idx);
ITR_daily = ITR(1:12);
ITR_cont = ITR(13:end-21);
X = [ones(length(ITR_daily),1),(session_break(perf_idx(1:length(ITR_daily))))'];
[b1,bint1,~,~,S1] = regress(ITR_daily,X);
X = [ones(length(ITR_cont),1),(session_break(perf_idx(13:end-21)))'];
[b2,bint2,~,~,S2] = regress(ITR_cont,X);


fprintf('\nPanel 2 Stats:\n')
fprintf('\nStat - ITR (first 5 days vs. last 5 days):\n')
fprintf('t-test - eff_sz=%.03f,   p=%.03f\n',mu1,p1)

% PLOTS - PANEL 2
subplot(4,1,2)
hold on
ylim([0.1,0.55])
% ylim([0,0.6])
xlim([0,trials(end)+1])
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color',[.7,.7,.7]);
end
scatter(session_break(perf_idx),perf{17}(idx),'o','markeredgecolor','w',...
    'sizedata',sz,'markerfacecolor','k'); hold on
plot(session_break(perf_idx),smooth(session_break(perf_idx),perf{17}(idx),20),'--k')
ylabel('bits / sec')
title('Fitt''s ITR')
set(gca,'XTick',[])

%% PLOTS: PANEL 3 & individual features

% look at angle btw velocity vectors (relative to final velocity within day)
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:));
Cy = squeeze(C(:,4,:));
Cc = squeeze(C(:,5,:));

xy_ang = [];
xy_diff_ang=[];

for i=1:length(day_break)-1,
    ref_idx = day_break(i+1)-1;
    idx = day_break(i):ref_idx-1;
    xy_ref_vectors=[Cx(:,ref_idx),Cy(:,ref_idx),zeros(896,1)];
    xy_ref_ang=atan2d(Cy(:,ref_idx),Cx(:,ref_idx));
    
    for t=idx,
        xy_ang= atan2d(Cy(:,t),Cx(:,t));
        xy_diff_ang(end+1) = rad2deg(subspace(xy_ang,xy_ref_ang));
        
    end
    xy_diff_ang(end+1) = nan;

end


% for all features together
figure('position',[445 326 1088 105]);
%set(gca,'ColorOrder',cc(4:6,:))
hold on
plot(xy_diff_ang,'linewidth',1);
% h(3)=plot(c_ang);
ylabel('Angle btw Decoders')
title('Within Day Decoder Convergence')
xlim([0,trials(end)+1])
ylim([0,90])
set(gca,'XTick',[],'YTick',[0,45,90])

for i=1:length(day_break)-2,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
for i=1:length(session_break)-1,
    vline(gca,session_break(i+1)-1,'color',[.7,.7,.7]);
end
% for i=1:length(days),
%     text(day_break(i),80,sprintf('%s/%s',...
%         days{i}(6),days{i}(7:8)))
% end
% for i=1:length(session_init),
%     if session_init(i)==3,
%         plot(session_break(i)+10,10,'*r')
%     end
% end
set(gca,'Children',flipud(get(gca,'Children')))


%% PLOTS: PANEL 4 & individual features
% look at angle btw velocity vectors (relative to final velocity across days)
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:));
Cy = squeeze(C(:,4,:));
Cc = squeeze(C(:,5,:));

% All features together
xy_ang = [];
xy_diff_ang=[];

ref_idx = trials(end);
idx = 1:ref_idx-1;
xy_ref_ang=atan2d(Cy(:,ref_idx),Cx(:,ref_idx));

for t=idx,
    xy_ang= atan2d(Cy(:,t),Cx(:,t));
    xy_diff_ang(end+1) = rad2deg(subspace(xy_ang,xy_ref_ang));
    
end
xy_diff_ang(end+1) = nan;

% All features together
figure('position',[445 326 1030 120]);
%set(gca,'ColorOrder',cc(4:6,:))
%hold on
plot(xy_diff_ang,'linewidth',1);
%h(2)=plot(y_ang,'linewidth',1);
% h(3)=plot(c_ang);
xlabel('trials')
title('Across Day Decoder Convergence')
xlim([0,trials(end)+1])
ylim([0,90])
set(gca,'YTick',[0,45,90])

for i=length(session_break)-1,
    vline(gca,session_break(i+1)-1,'color','k','linewidth',2);
end
for i=1:length(session_break)-1,
    vline(gca,session_break(i+1)-1,'color',[.7,.7,.7]);
end

set(gca,'Children',flipud(get(gca,'Children')))
% legend(h,{'Vx','Vy','C'},'location','northeast')
%legend(h,{'Vx','Vy'},'location','northeast')


% for each feature: calculation
xy_ang = [];
xy_diff_ang=[];
xy_diff_ang_perFe=[];

ref_idx = trials(end);
idx = 1:ref_idx-1;
xy_ref_ang=atan2d(Cy(:,ref_idx),Cx(:,ref_idx));

for Fe=1:7
    Index_Fe=128*(Fe-1)+1:(128*(Fe-1)+128);
    Counter=0;
    for t=idx,
        Counter=Counter+1;
        xy_ang= atan2d(Cy(Index_Fe,t),Cx(Index_Fe,t));
        xy_diff_ang_perFe(Counter,Fe) = rad2deg(subspace(xy_ang,xy_ref_ang(Index_Fe,1)));
        
    end
end

% plotting 
FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};
figure('position',[445 326 1300 400]);
for Fe=1:7
    plot(xy_diff_ang_perFe(:,Fe),'linewidth',2);
    %legend(FeatureName{Fe});
    hold on
    
end 
xlabel('trials')
title('Across Day Vectorial Convergence of Decoder Per Feature')
xlim([0,trials(end)+1])
ylim([0,100])
set(gca,'YTick',[0,45,90])
legend(FeatureName{1},FeatureName{2},FeatureName{3},FeatureName{4},FeatureName{5},...
    FeatureName{6},FeatureName{7});
%HighQualityFigs('Convergence_Vectorial_PerFeature')






