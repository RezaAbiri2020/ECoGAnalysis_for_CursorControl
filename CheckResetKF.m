%% Investigate the variation in KF weights after reset using imagined blocks with three features set delta, beta, gamma2
% the CLDA was performed in Center-Out task
% calculate variation in KF weights 

clear all 
close all  


%% Experiment Info 
% for the following day: CursorControlGridTask
if 0
  expts(end+1).yymmdd = '20200127';
expts(end).CLDA_hhmmss = {'114408','140423'};
expts(end).perf_CenterOut_hhmmss = {'112639','140423'};
expts(end).perf_RadialTask_hhmmss = {'141809','142333','142707'};  
end 

if 0
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20200122';
exptsCLDA(end).hhmmss = {'113333','145510','151403'};
 

% for the following day: CursorControlGridTask
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20200124';
exptsCLDA(end).hhmmss = {'140515'};
 
% for the following day: CursorControlGridTask
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20200127';
exptsCLDA(end).hhmmss = {'114408','140423'};


exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20200131';
exptsCLDA(end).hhmmss = {'113822','133056'};
end 

exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20200204';
exptsCLDA(end).hhmmss = {'114712','133740'};




%% Plotting the convergence of reset KF compared to PnP 20191217 for all 3 features and individual ones

% going through CLDA files and find the angle with reference
Angle_x_all = [];
Angle_x_PerFe=[];
Angle_y_all = [];
Angle_y_PerFe=[]; 
Angle_c_all = [];
Angle_c_PerFe=[];
Sessions={};

% the ref
datadir = fullfile('Z:\Bravo\Bravo1\','20191217',...
    'GangulyServer','20191217','CenterOut','134525','BCI_CLDA');
files = dir(fullfile(datadir,'Data*.mat'));
load(fullfile(datadir,files(end).name));
C_Ref=TrialData.KalmanFilter{1,1}.C;

for i=1:length(exptsCLDA)
    expt = exptsCLDA(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.hhmmss)
        
        % go through datafiles in CLDA blocks
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer',yymmdd,'CenterOut',expt.hhmmss{1,j},'BCI_CLDA');
            
       
        Sessions(end+1)={[yymmdd,',',expt.hhmmss{1,j}]};
            
        files = dir(fullfile(datadir,'Data*.mat'));
        
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            C_New=TrialData.KalmanFilter{1,1}.C;
            Angle_x_all = [Angle_x_all;subspace(C_New(:,3),C_Ref(:,3))*180/pi];
            Angle_y_all = [Angle_y_all;subspace(C_New(:,4),C_Ref(:,4))*180/pi];
            Angle_c_all = [Angle_c_all;subspace(C_New(:,5),C_Ref(:,5))*180/pi];
            
            for feature=1:3
                AnglesX(feature)= subspace(C_New(128*(feature-1)+1:128*feature,3),C_Ref(128*(feature-1)+1:128*feature,3))*180/pi;
                AnglesY(feature) = subspace(C_New(128*(feature-1)+1:128*feature,4),C_Ref(128*(feature-1)+1:128*feature,4))*180/pi;
                AnglesC(feature) = subspace(C_New(128*(feature-1)+1:128*feature,5),C_Ref(128*(feature-1)+1:128*feature,5))*180/pi;   
            end
            Angle_x_PerFe=[Angle_x_PerFe; AnglesX];
            Angle_y_PerFe=[Angle_y_PerFe; AnglesY];
            Angle_c_PerFe=[Angle_c_PerFe; AnglesC];
            
            
        end
        
        clear files
    end
    
    
end
% for all fearures
figure('position',[400 400 1800 800]);
plot(Angle_x_all,'linewidth',2)
hold on
plot(Angle_y_all,'linewidth',2)
hold on
plot(Angle_c_all,'linewidth',2)
set(gca,'FontSize',20)
legend ('Vx','Vy','Const','Location','Northeast')
ylabel('Degree','FontSize',14)
xlabel('CLDA Trials','FontSize',14)
title(['the process of converging ResetKF (from imagined Mov:',exptsCLDA.yymmdd,') to PnP decoder:20191217'])
HighQualityFigs(['ResetKF_',exptsCLDA.yymmdd,'_AllFeatures'])

% for Vx per feature
figure('position',[400 400 1800 800]);
plot(Angle_x_PerFe,'linewidth',2)
set(gca,'FontSize',20)
legend ('delta','beta','hg','Location','Northeast')
ylabel('Degree','FontSize',14)
xlabel('CLDA Trials','FontSize',14)
title(['Vx:the process of converging ResetKF (from imagined Mov:',exptsCLDA.yymmdd,') to PnP decoder:20191217'])
HighQualityFigs(['ResetKF_Vx_',exptsCLDA.yymmdd,'_PerFeature'])

% for Vy per feature
figure('position',[400 400 1800 800]);
plot(Angle_y_PerFe,'linewidth',2)
set(gca,'FontSize',20)
legend ('delta','beta','hg','Location','Northeast')
ylabel('Degree','FontSize',14)
xlabel('CLDA Trials','FontSize',14)
title(['Vy:the process of converging ResetKF (from imagined Mov:',exptsCLDA.yymmdd,') to PnP decoder:20191217'])
HighQualityFigs(['ResetKF_Vy_',exptsCLDA.yymmdd,'_PerFeature'])

% for constant term per feature
figure('position',[400 400 1800 800]);
plot(Angle_c_PerFe,'linewidth',2)
set(gca,'FontSize',20)
legend ('delta','beta','hg','Location','Northeast')
ylabel('Degree','FontSize',14)
xlabel('CLDA Trials','FontSize',14)
title(['Const:the process of converging ResetKF (from imagined Mov:',exptsCLDA.yymmdd,') to PnP decoder:20191217'])
HighQualityFigs(['ResetKF_C_',exptsCLDA.yymmdd,'_PerFeature'])



%% plot the trend of changes across CLDA trials for feature weights using C martix in maps
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

% the ref map
datadir = fullfile('Z:\Bravo\Bravo1\','20191217',...
    'GangulyServer','20191217','CenterOut','134525','BCI_CLDA');
files = dir(fullfile(datadir,'Data*.mat'));
load(fullfile(datadir,files(end).name));
C_Ref=TrialData.KalmanFilter{1,1}.C;

figure('position',[400 400 1800 300]);
%suptitle('The magnitude of features weights: PnP: 20191217')
for feature=1:3
    Cx_feature = C_Ref(128*(feature-1)+1:128*feature,3);
    Cy_feature = C_Ref(128*(feature-1)+1:128*feature,4);
    cc = brewermap(75,'YlOrRd');
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

%HighQualityFigs('Maps_PnP20191217')

Sessions={};

for i=1:length(exptsCLDA)
    expt = exptsCLDA(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.hhmmss)
        
        % go through datafiles in CLDA blocks
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer',yymmdd,'CenterOut',expt.hhmmss{1,j},'BCI_CLDA');
            
        Sessions(end+1)={[yymmdd,'-',expt.hhmmss{1,j}]};
            
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
            HighQualityFigs(['Maps_',yymmdd,'_CLDA',num2str(j),'-',num2str(kk)])
        end
        
        
        
    end 
    
end 

%%
