
%%
% to address the comments related to PnP. In order to adress this comments,
% new terminology is defined such as Micro-CLDA

%% Investigate the variation in KF weights for Micro-CLDA in Center-Out task across multiple months with full-feature set
clear all 
close all

% experiment info 
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20190724';
exptsCLDA(end).hhmmss = {'110736','140225'};

exptsCLDA(end+1).yymmdd = '20190725';
exptsCLDA(end).hhmmss = {'112320','135821'};

exptsCLDA(end+1).yymmdd = '20190726';
exptsCLDA(end).hhmmss = {'112855','131757'};

exptsCLDA(end+1).yymmdd = '20190730';
exptsCLDA(end).hhmmss = {'112918','140229','145137'};

exptsCLDA(end+1).yymmdd = '20190807';
exptsCLDA(end).hhmmss = {'111451','135559'};

exptsCLDA(end+1).yymmdd = '20190809';
exptsCLDA(end).hhmmss = {'112404','140130'};

exptsCLDA(end+1).yymmdd = '20190816';
exptsCLDA(end).hhmmss = {'110820'};

exptsCLDA(end+1).yymmdd = '20190830';
exptsCLDA(end).hhmmss = {'112710','140606'};

exptsCLDA(end+1).yymmdd = '20190904';
exptsCLDA(end).hhmmss = {'113405','140826'};

exptsCLDA(end+1).yymmdd = '20190906';
exptsCLDA(end).hhmmss = {'110649','134415'};

% going through CLDA files and find the angle with reference
Angle_x = [];
Angle_y = [];
Angle_c = [];
Sessions={};

for i=1:length(exptsCLDA)
    expt = exptsCLDA(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.hhmmss)
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.hhmmss{1,j},'BCI_CLDA');
        Sessions(end+1)={[yymmdd,',',expt.hhmmss{1,j}]};
            
        files = dir(fullfile(datadir,'Data*.mat'));
        load(fullfile(datadir,files(end).name));
        
        if i==1 && j==1
            C_Ref=TrialData.KalmanFilter{1,1}.C;
        end
        
        C_New=TrialData.KalmanFilter{1,1}.C;
        
        Angle_x = [Angle_x;subspace(C_New(:,3),C_Ref(:,3))*180/pi];
        Angle_y = [Angle_y;subspace(C_New(:,4),C_Ref(:,4))*180/pi];
        Angle_c = [Angle_c;subspace(C_New(:,5),C_Ref(:,5))*180/pi];
        
    clear files
    end
    
    
end

figure('position',[400 400 1800 800]);
plot(Angle_x,'linewidth',2)
hold on
plot(Angle_y,'linewidth',2)
hold on
plot(Angle_c,'linewidth',2)
set(gca,'FontSize',20)
legend ('Vx','Vy','Const','Location','Northwest')
ylabel('Degree','FontSize',14)
xticks([1:1:length(Angle_x)])
xtickangle(45)
xticklabels(Sessions)
HighQualityFigs('MicroCLDA')


%% Investigate the variation in KF weights for PnP in radial task with three features set delta, beta, gamma2
% the CLDA was performed in Center-Out task
% calculate variation in KF weights 

clear all 
close all 

% experiment info for the first fixed blocks of PnP 

% for the following day: CursorControlGridTask
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20191119';
exptsCLDA(end).hhmmss = {'111123','111851'};

% for the following day: CursorControlGridTask
exptsCLDA(end+1).yymmdd = '20191125';
exptsCLDA(end).hhmmss = {'113919','114901','151054'};

% for the following day: RadialTypingMultiClick
exptsCLDA(end+1).yymmdd = '20191127';
exptsCLDA(end).hhmmss = {'135842','140144'};

% for the following day: RadialTask
exptsCLDA(end+1).yymmdd = '20200115';
exptsCLDA(end).hhmmss = {'110726','111243'};

% for the following day: RadialTask
exptsCLDA(end+1).yymmdd = '20200117';
exptsCLDA(end).hhmmss = {'112153','113322','114206'};


% going through CLDA files and find the angle with reference
Angle_x = [];
Angle_y = [];
Angle_c = [];
Sessions={};

for i=1:length(exptsCLDA)
    expt = exptsCLDA(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.hhmmss)
        
        % go through datafiles in CLDA blocks
        if i<=2
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer','CursorControlGridTask',yymmdd,expt.hhmmss{1,j},'BCI_Fixed');
        elseif i==3
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'RadialTypingMultiClick',expt.hhmmss{1,j},'BCI_Fixed');
        else
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer',yymmdd,'RadialTask',expt.hhmmss{1,j},'BCI_Fixed');
            
        end
        Sessions(end+1)={[yymmdd,',',expt.hhmmss{1,j}]};
            
        files = dir(fullfile(datadir,'Data*.mat'));
        load(fullfile(datadir,files(end).name));
        
        if i==1 && j==1
            C_Ref=TrialData.KalmanFilter{1,1}.C;
        end
        
        C_New=TrialData.KalmanFilter{1,1}.C;
        
        Angle_x = [Angle_x;subspace(C_New(:,3),C_Ref(:,3))*180/pi];
        Angle_y = [Angle_y;subspace(C_New(:,4),C_Ref(:,4))*180/pi];
        Angle_c = [Angle_c;subspace(C_New(:,5),C_Ref(:,5))*180/pi];
        
    clear files
    end
    
    
end

figure('position',[400 400 1800 800]);
plot(Angle_x,'linewidth',2)
hold on
plot(Angle_y,'linewidth',2)
hold on
plot(Angle_c,'linewidth',2)
set(gca,'FontSize',20)
legend ('Vx','Vy','Const','Location','Northwest')
ylabel('Degree','FontSize',14)
xticks([1:1:length(Angle_x)])
xtickangle(45)
xticklabels(Sessions)
HighQualityFigs('RadialFixed')



%% checkup for which KF weights we used for plug and play 
clear all 
close all 

% experiment info

% for the following day: CenterOut
exptsCLDA = [];
exptsCLDA(end+1).yymmdd = '20191217';
exptsCLDA(end).hhmmss = {'111414','134525'};

% for the following day: RadialTask
exptsCLDA(end+1).yymmdd = '20200115';
exptsCLDA(end).hhmmss = {'110726','111243'};

% for the following day: RadialTask
exptsCLDA(end+1).yymmdd = '20200117';
exptsCLDA(end).hhmmss = {'112153','113322','114206'};


% going through CLDA files and find the angle with reference
Angle_x = [];
Angle_y = [];
Angle_c = [];
Sessions={};

for i=1:length(exptsCLDA)
    expt = exptsCLDA(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.hhmmss)
        
        % go through datafiles in CLDA blocks
        if i==1
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer',yymmdd,'CenterOut',expt.hhmmss{1,j},'BCI_Fixed');
        else
            datadir = fullfile('Z:\Bravo\Bravo1\',yymmdd,...
                'GangulyServer',yymmdd,'RadialTask',expt.hhmmss{1,j},'BCI_Fixed');
            
        end
        Sessions(end+1)={[yymmdd,',',expt.hhmmss{1,j}]};
            
        files = dir(fullfile(datadir,'Data*.mat'));
        load(fullfile(datadir,files(end).name));
        
        if i==1 && j==1
            C_Ref=TrialData.KalmanFilter{1,1}.C;
        end
        
        C_New=TrialData.KalmanFilter{1,1}.C;
        
        Angle_x = [Angle_x;subspace(C_New(:,3),C_Ref(:,3))*180/pi];
        Angle_y = [Angle_y;subspace(C_New(:,4),C_Ref(:,4))*180/pi];
        Angle_c = [Angle_c;subspace(C_New(:,5),C_Ref(:,5))*180/pi];
        
    clear files
    end
    
    
end

figure('position',[400 400 1800 800]);
plot(Angle_x,'linewidth',2)
hold on
plot(Angle_y,'linewidth',2)
hold on
plot(Angle_c,'linewidth',2)
set(gca,'FontSize',20)
legend ('Vx','Vy','Const','Location','Northwest')
ylabel('Degree','FontSize',14)
xticks([1:1:length(Angle_x)])
xtickangle(45)
xticklabels(Sessions)
HighQualityFigs('KFWeights_ForPnP')


