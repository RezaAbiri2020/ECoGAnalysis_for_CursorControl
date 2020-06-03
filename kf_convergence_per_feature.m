%% check if kalman weights are converging within and across days, compare
% with performance
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

%% combine all plots into single figure

% look at angle btw weights (relative to final weights across days)
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:));
Cy = squeeze(C(:,4,:));
Cc = squeeze(C(:,5,:));


ref_idx = trials(end);
idx = 1:ref_idx-1;

figure('position',[445 326 1088 547]);
subplot(4,1,4); hold on
% cc1 = brewermap(7,'*BuGn');
% cc2 = brewermap(7,'GnBu');
% cc = cat(1,cc1(1:3,:),cc2(4:end,:));
cc = brewermap(3,'*Dark2');
% cc = get(groot,'defaultAxesColorOrder');

fs = [2,4,7];
for f=1:3,%1:7,
    f_idx = (1:128) + 128*(fs(f)-1); % feature idx
    
    x_ang = [];
    y_ang = [];
    c_ang = [];
    for t=idx,
        x_ang(end+1) = rad2deg(subspace(Cx(f_idx,t),Cx(f_idx,ref_idx)));
        y_ang(end+1) = rad2deg(subspace(Cy(f_idx,t),Cy(f_idx,ref_idx)));
        c_ang(end+1) = rad2deg(subspace(Cc(f_idx,t),Cc(f_idx,ref_idx)));
    end
    x_ang(end+1) = nan;
    y_ang(end+1) = nan;
    c_ang(end+1) = nan;
    
    % add line to plot
    h(f)=plot(.5*(x_ang+y_ang),'linewidth',1,'color',cc(f,:));
    
end

% h(3)=plot(c_ang);
xlabel('trials')
title('Across Day Decoder Convergence')
xlim([0,trials(end)+1])
ylim([0,90])
set(gca,'YTick',[0,45,90])

for i=length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color','k','linewidth',2);
end
for i=1:length(day_break)-1,
    vline(gca,day_break(i+1)-1,'color',[.7,.7,.7]);
end

set(gca,'Children',flipud(get(gca,'Children')))
% legend(h,{'\angle\delta','|\delta|','|\theta|','|\alpha|',...
%     '|\beta|','|\gamma_1|','|\gamma_2|'},'location','northeast')
legend(h,{'|\delta|','|\mu|','|\gamma_2|'},'location','northeast')

% save figure
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/kalman_C_angle_across_days_per_feature.pdf',...
%     '-pdf','-r300','-painters')



%% stats: angle for each feature, 
% regression of angle in first 5 days and middle 5 days
f_ang_all = [];
for f=1:7,
    f_idx = (1:128) + 128*(f-1); % feature idx
    
    x_ang = [];
    y_ang = [];
    for t=idx,
        x_ang(end+1) = rad2deg(subspace(Cx(f_idx,t),Cx(f_idx,ref_idx)));
        y_ang(end+1) = rad2deg(subspace(Cy(f_idx,t),Cy(f_idx,ref_idx)));
    end
    f_ang_all(:,f) = .5*(x_ang(:) + y_ang(:));
end

idx_daily = 1:day_break(5); % 1st 5 days are daily
f_ang_daily = f_ang_all(idx_daily,:);

idx_cont = day_break(5)+1:day_break(10); % middle 5 days are cont
f_ang_cont = f_ang_all(idx_cont,:);

fs = [1,2,3,4,5,6,7];
b_daily = {}; bint_daily = {}; S_daily = {};
b_cont = {}; bint_cont = {}; S_cont = {};
for i=1:7,
    f = fs(i);
    N = size(f_ang_daily,1);
    X = [ones(N,1),(1:N)'];
    [b_daily{i},bint_daily{i},~,~,S_daily{i}] = regress(f_ang_daily(:,f),X,0.05/3);

    N = size(f_ang_cont,1);
    X = [ones(N,1),(1:N)'];
    [b_cont{i},bint_cont{i},~,~,S_cont{i}] = regress(f_ang_cont(:,f),X,0.05/3);
end

% print stats
fprintf('\nSTATS: regression of feature decoder angles wrt final decoder\n')
fprintf('  reporting regressions of continuous (1st 5 days)\n')
fprintf('   feature convergence rates\n')
for i=1:7,
    f = fs(i);
    fprintf('  \nFeature%i: int=%.03f, slope=%.03f, CI=[%.03f,%.03f]\n',...
        f,b_cont{i}(1),b_cont{i}(2),bint_cont{i}(2,1),bint_cont{i}(2,2))
end

