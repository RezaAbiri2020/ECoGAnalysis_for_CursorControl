%% check if kalman weights are converging within and across days
clear, clc, close all

%% experiment info
expts = [];

expts(end+1).yymmdd = '20190501';
expts(end).hhmmss = {'133745'};

expts(end+1).yymmdd = '20190507';
expts(end).hhmmss = {'112817'};

expts(end+1).yymmdd = '20190510';
expts(end).hhmmss = {'140749'};

expts(end+1).yymmdd = '20190514';
expts(end).hhmmss = {'133050'};

expts(end+1).yymmdd = '20190515';
expts(end).hhmmss = {'110735'};

expts(end+1).yymmdd = '20190521';
expts(end).hhmmss = {'141438'};

expts(end+1).yymmdd = '20190524';
expts(end).hhmmss = {'133313'};

expts(end+1).yymmdd = '20190529';
expts(end).hhmmss = {'113247'};

expts(end+1).yymmdd = '20190531';
expts(end).hhmmss = {'135046'};

expts(end+1).yymmdd = '20190604';
expts(end).hhmmss = {'140706'};

expts(end+1).yymmdd = '20190607';
expts(end).hhmmss = {'140114'};

expts(end+1).yymmdd = '20190610';
expts(end).hhmmss = {'104104'};

expts(end+1).yymmdd = '20190618';
expts(end).hhmmss = {'141917'};

expts(end+1).yymmdd = '20190621';
expts(end).hhmmss = {'135917'};

expts(end+1).yymmdd = '20190626';
expts(end).hhmmss = {'112324'};


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

        for iii=T,
            % load data, grab neural features
            disp(datafiles(iii).name)
            load(fullfile(datadir,datafiles(iii).name))
            
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

%% figure info 
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
limch = ch_layout(R,1);

% go through each feature and plot erps
feature_strs = {'delta-phase','delta-pwr','theta-pwr','mu-pwr',...
    'beta-pwr','low-gamma-pwr','high-gamma-pwr'};
def_cc = get(groot,'defaultAxesColorOrder');
cc = brewermap(trials(end),'Blues');
cc(1:5,:) = repmat(def_cc(2,:),5,1);

%% plot kalman preferred dirs
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';

% features
for feature=1:7,
    Cx_feature = Cx(:,128*(feature-1)+1:128*feature);
    Cy_feature = Cy(:,128*(feature-1)+1:128*feature);

    fig = figure('units','normalized',...
        'position',[0.1010 0.5102 0.3281 0.3380],...
        'name',feature_strs{feature});
    set(gca,'NextPlot','add');
    xx = [min(Cx_feature(end,:)),max(Cx_feature(end,:))];
    yy = [min(Cy_feature(end,:)),max(Cy_feature(end,:))];
    xy = [min([xx,yy]),max([xx,yy])];
    w = mean(abs(xy));
%     w = 6e-4; % manual for low gamma fig
%     w = 25e-4; % manual for high gamma fig
    for sess=1:trials(end),
        Cvec = [Cx(sess,128*(feature-1)+1:128*feature);
            Cy(sess,128*(feature-1)+1:128*feature)]';

        for ch=1:Nch,
            [r,c] = find(ch_layout == ch);
            x = w*c;
            y = w*(R-r);

            % clean plot
            plot(x+[0,Cvec(ch,1)],y+[0,Cvec(ch,2)],...
                'linewidth',0.5,'color',cc(sess,:))
        end
    end

    title(strrep(...
        sprintf('kalman_%s_pref_dir_across_days',...
        feature_strs{feature}),...
        '_',' '))
    axis equal
    box on
    
    % save figure
%     export_fig(sprintf('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/KF_PD/kalman_%s_pref_dir_across_days.pdf',...
%         feature_strs{feature}),'-pdf','-r300','-painters')
%     close(fig)

end

%% plot kalman preferred dirs - just a few example channels on cartesian plots
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';

% features
for feature=1:7,
    Cx_feature = Cx(:,128*(feature-1)+1:128*feature);
    Cy_feature = Cy(:,128*(feature-1)+1:128*feature);

    xx = [min(Cx_feature(end,:)),max(Cx_feature(end,:))];
    yy = [min(Cy_feature(end,:)),max(Cy_feature(end,:))];
    xy = [min([xx,yy]),max([xx,yy])];
    w = mean(abs(xy));
%     w = 6e-4; % manual for low gamma fig
%     w = 25e-4; % manual for high gamma fig

    fig = figure('units','normalized',...
        'position',[0.0823 0.3935 0.4615 0.4935],...
        'name',feature_strs{feature});
    chs = [114,122,80,52,11];
    for i=1:length(chs),
        ch = chs(i);
        subplot(2,3,i)
        set(gca,'NextPlot','add');
        
        for sess=1:trials(end),
            Cvec = [Cx(sess,128*(feature-1)+1:128*feature);
                Cy(sess,128*(feature-1)+1:128*feature)]';

            x = [0,Cvec(ch,1)];
            y = [0,Cvec(ch,2)];
            plot(x,y,'color',cc(sess,:))
            
        end
        axis equal
        set(gca,'xlim',[-15e-4,2e-4],'ylim',[-2e-4,18e-4])
%         set(gca,'xlim',[-5e-4,5e-4],'ylim',[-5e-4,10e-4])
        title(sprintf('ch%i',ch))
        
    end
    
    % save figure
%     export_fig(sprintf('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/KF_PD/kalman_%s_pref_dir_across_days_examples_v2.pdf',...
%         feature_strs{feature}),'-pdf','-r300','-painters')

end

%% plot norm of each feature weight as an image
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';
cc = brewermap(75,'YlOrRd');

% features
for feature=6%:7,
    figure('units','normalized','position',[0.3516 0.6065 0.3630 0.2944]);
    
    for sess=1:trials(end),
        Cx_feature = Cx(:,128*(feature-1)+1:128*feature);
        Cy_feature = Cy(:,128*(feature-1)+1:128*feature);
        
        Cvec = [Cx(sess,128*(feature-1)+1:128*feature);
            Cy(sess,128*(feature-1)+1:128*feature)]';
        
        % magnitude of kalman weights
        mag = zeros(size(Cvec,1),1);
        for i=1:size(Cvec,1),
            mag(i) = norm(Cvec(i,:));
        end
        
        % put into plotting matrix
        x = zeros(size(ch_layout));
        y = zeros(size(ch_layout));
        M = zeros(size(ch_layout));
        for ch=1:Nch,
            [r,c] = find(ch_layout == ch);
            x(r,c) = c;
            y(r,c) = (R-r)+1;
            M(r,c) = mag(ch);
        end
        
        % plot
        ax(sess)=subplot(3,5,sess);
        imagesc(M)
        colormap(cc)
        title(sprintf('%s',days{sess}))
        
    end % day
    
    % clean plot
%     z = cell2mat(get(ax,'CLim'));
%      zz = [min(z(:)),max(z(:))];

%      zz = [0,0.0007]; % manual for beta
%    zz = [0,0.0005]; % manual for mu
     zz = [0,0.00075]; % manual for low gamma
%     zz = [0,0.003]; % manual for high gamma

    set(ax,'CLim',zz,'XTick',[],'YTick',[])
    suptitle(sprintf('Kalman Weight Magnitude: %s',...
        feature_strs{feature}))
    
    % add colorbar to last plot
    colorbar(ax(8),'location','southoutside','position',[.4,.06,.2,.025])
    
    % save
%     export_fig(sprintf('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/KF_PD_Mag/kalman_%s_C_mag_across_days_image.pdf',...
%         feature_strs{feature}),'-pdf','-r300','-painters')
%     close(gcf)
    
end % feature

%% plot norm of each feature weight as an image on final day to compare feature maps
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';
cc = brewermap(75,'YlOrRd');

figure('units','normalized','position',[0.1328 0.5222 0.6214 0.2519]);

% features
for feature=1:7,
    
    for sess=trials(end),
        Cx_feature = Cx(:,128*(feature-1)+1:128*feature);
        Cy_feature = Cy(:,128*(feature-1)+1:128*feature);
        
        Cvec = [Cx(sess,128*(feature-1)+1:128*feature);
            Cy(sess,128*(feature-1)+1:128*feature)]';
        
        % magnitude of kalman weights
        mag = zeros(size(Cvec,1),1);
        for i=1:size(Cvec,1),
            mag(i) = norm(Cvec(i,:));
        end
        
        % put into plotting matrix
        x = zeros(size(ch_layout));
        y = zeros(size(ch_layout));
        M = zeros(size(ch_layout));
        for ch=1:Nch,
            [r,c] = find(ch_layout == ch);
            x(r,c) = c;
            y(r,c) = (R-r)+1;
            M(r,c) = mag(ch);
        end
        
        % plot
        ax=subplot(2,4,feature);
        imagesc(M)
        colormap(cc)
        title(sprintf('%s',feature_strs{feature}))
        
    end % day
    
    % clean plot
%     z = cell2mat(get(ax,'CLim'));
%     zz = [min(z(:)),max(z(:))];
%     zz = [0,0.00075]; % manual for low gamma
%     zz = [0,0.003]; % manual for high gamma
    set(ax,'XTick',[],'YTick',[])
%     suptitle('Kalman Weight Maps')
    
    % add colorbar to last plot
    colorbar(ax,'location','westoutside')
    
    
    
end % feature
% save
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/KF_PD_Mag/kalman_all_features_C_mag_image.pdf',...
%     '-pdf','-r300','-painters')
% close(gcf)

%% plot trend of each feature weight mag across chs - hg
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';
def_cc = get(groot,'defaultAxesColorOrder');
cc = brewermap(trials(end),'Blues');
cc(1:5,:) = repmat(def_cc(2,:),5,1);
alpha = [.5*ones(1,5),ones(1,10)];
edges = 0:1e-4:25e-4;
binc = (edges(1:end-1) + edges(2:end))/2;

% features
M = [];
for feature=7,
    
    for sess=1:trials(end),
        
        % magnitude of kalman weights
        Cvec = [Cx(sess,:);Cy(sess,:)]';
        mag = zeros(size(Cvec,1),1);
        for i=1:size(Cvec,1),
            mag(i) = norm(Cvec(i,:));
        end
        
        
        % grab all chans for feature and discretize
        M(:,sess) = mag(128*(feature-1)+1:128*feature);
        
    end % day
    
    % get slopes
    b = M/[ones(1,15);(1:15)];
    slopes = b(:,2);
    
    % plot distribution
    figure('position',[665 774 348 163])
    ax = gca;
    hold on
    histogram(slopes,-10.5e-5:1e-5:10.5e-5,'normalization','probability',...
        'orientation','vertical')
    ylabel('slope across channels')
    xlabel('pdf')
    
    neg_ch = 128;
    pos_ch = 114;
    vline(gca,b(neg_ch,2),'-');
    vline(gca,b(pos_ch,2),'-');
    
    % plot example trend line (positive)
    ch = neg_ch; % 106, 114
    axes('position',[.2,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\gamma_2|_{ch%i}',ch))
    

    % plot example trend line (negative)
    ch = pos_ch; % 96, 128
    axes('position',[.65,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\gamma_2|_{ch%i}',ch))
    
    % report stats for pos ch
    [b,bint] = regress(M(pos_ch,:)',[ones(15,1),(1:15)']);
    fprintf('\nSTATS: lin regression (gamma2 pos ch):\n')
    fprintf('  int=%.03g, slope=%.03g, CI=[%.03g,%.03g]\n',...
        b(1),b(2),bint(2,1),bint(2,2))
    
end % feature

% save
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/kalman_hg_mag_trends_across_days_distr.pdf',...
%     '-pdf','-r300','-painters')
% close(gcf)

%% plot trend of each feature weight mag across chs - mu
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';
def_cc = get(groot,'defaultAxesColorOrder');
cc = brewermap(trials(end),'Blues');
cc(1:5,:) = repmat(def_cc(2,:),5,1);
alpha = [.5*ones(1,5),ones(1,10)];
edges = 0:1e-4:25e-4;
binc = (edges(1:end-1) + edges(2:end))/2;

% features
for feature=4,
    
    for sess=1:trials(end),
        
        % magnitude of kalman weights
        Cvec = [Cx(sess,:);Cy(sess,:)]';
        mag = zeros(size(Cvec,1),1);
        for i=1:size(Cvec,1),
            mag(i) = norm(Cvec(i,:));
        end
        
        
        % grab all chans for feature and discretize
        M(:,sess) = mag(128*(feature-1)+1:128*feature);
        
    end % day
    
    % get slopes
    b = M/[ones(1,15);(1:15)];
    slopes = b(:,2);
    
    % plot distribution
    figure('position',[665 774 348 163])
    ax = gca;
    hold on
    histogram(slopes,-5e-5:.5e-5:2e-5,'normalization','probability',...
        'orientation','vertical')
    ylabel('slope across channels')
    xlabel('pdf')
    
    neg_ch = 53;
    pos_ch = 26;
    vline(gca,b(neg_ch,2),'-');
    vline(gca,b(pos_ch,2),'-');
    
    % plot example trend line (positive)
    ch = neg_ch; % 106, 114
    axes('position',[.2,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\mu|_{ch%i}',ch))
    

    % plot example trend line (negative)
    ch = pos_ch; % 96, 128
    axes('position',[.65,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\mu|_{ch%i}',ch))
    
    % report stats for neg ch
    [b,bint] = regress(M(neg_ch,:)',[ones(15,1),(1:15)']);
    fprintf('\nSTATS: lin regression (mu neg ch):\n')
    fprintf('  int=%.03g, slope=%.03g, CI=[%.03g,%.03g]\n',...
        b(1),b(2),bint(2,1),bint(2,2))
    
end % feature

% save
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/kalman_mu_mag_trends_across_days_distr.pdf',...
%     '-pdf','-r300','-painters')
% close(gcf)

%% plot trend of each feature weight mag across chs - delta
C = cat(3,Ccell{:});
Cx = squeeze(C(:,3,:))';
Cy = -1*squeeze(C(:,4,:))';
def_cc = get(groot,'defaultAxesColorOrder');
cc = brewermap(trials(end),'Blues');
cc(1:5,:) = repmat(def_cc(2,:),5,1);
alpha = [.5*ones(1,5),ones(1,10)];
edges = 0:1e-4:25e-4;
binc = (edges(1:end-1) + edges(2:end))/2;

% features
for feature=2,
    
    for sess=1:trials(end),
        
        % magnitude of kalman weights
        Cvec = [Cx(sess,:);Cy(sess,:)]';
        mag = zeros(size(Cvec,1),1);
        for i=1:size(Cvec,1),
            mag(i) = norm(Cvec(i,:));
        end
        
        
        % grab all chans for feature and discretize
        M(:,sess) = mag(128*(feature-1)+1:128*feature);
        
    end % day
    
    % get slopes
    b = M/[ones(1,15);(1:15)];
    slopes = b(:,2);
    
    % plot distribution
    figure('position',[665 774 348 163])
    ax = gca;
    hold on
    histogram(slopes,-10.5e-5:1e-5:10.5e-5,'normalization','probability',...
        'orientation','vertical')
    ylabel('slope across channels')
    xlabel('pdf')
    
    % delta
    neg_ch = 93;
    pos_ch = 16;
    vline(gca,b(neg_ch,2),'-');
    vline(gca,b(pos_ch,2),'-');
    
    % plot example trend line (positive)
    ch = neg_ch; % 106, 114
    axes('position',[.2,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\delta|_{ch%i}',ch))
    

    % plot example trend line (negative)
    ch = pos_ch; % 96, 128
    axes('position',[.65,.65,.25,.15]); hold on
    scatter(1:15,M(ch,:),20,cc,'o','filled','markeredgecolor','w')
    plot(1:15,b(ch,1) + b(ch,2)*(1:15), '-k')
    set(gca,'xtick',[])
    xlim([.5,15.5])
    xlabel('days')
    title(sprintf('|\\delta|_{ch%i}',ch))
    
end % feature

% save
% export_fig('~/Desktop/BCI_Learning_Manuscript/Figures/Decoder_Convergence/kalman_delta_mag_trends_across_days_distr.pdf',...
%     '-pdf','-r300','-painters')
% close(gcf)

