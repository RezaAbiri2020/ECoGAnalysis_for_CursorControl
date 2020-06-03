%% check if kalman weigths are converging across days and months
% also calculate the corresponding performance from fixed blocks 
% immidiately done after CLDA
clear all
close all
clc

%% experiment info
expts = [];

% Long-Term CLDA in Center-Out and test performance in Center-Out 
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

expts(end+1).yymmdd = '20190628';
expts(end).hhmmss = {'112947','113528'};
expts(end).perf_hhmmss = {'112947','113528','114111'};

% start of CursorControlGridTask
% Long-Term CLDA in Center-Out and test performance in square grid 
expts(end+1).yymmdd = '20190710';
expts(end).hhmmss = {'113803','134721'};
expts(end).perf_hhmmss = {'115126','115539','135416'};

expts(end+1).yymmdd = '20190712';
expts(end).hhmmss = {'133612'};
expts(end).perf_hhmmss = {'134505'};

expts(end+1).yymmdd = '20190726';
expts(end).hhmmss = {'131757'};
expts(end).perf_hhmmss = {'133853'};

% expts(end+1).yymmdd = '20190730';
% expts(end).hhmmss = {'140229','145137'};
% expts(end).perf_hhmmss = {'143832','145725'};

expts(end+1).yymmdd = '20190807';
expts(end).hhmmss = {'111451','135559'};
expts(end).perf_hhmmss = {'112500','113225','140450','141204'};

expts(end+1).yymmdd = '20190809';
expts(end).hhmmss = {'112404','140130'};
%, excluding this run: '141115'
expts(end).perf_hhmmss = {'114432','114630','114836','142846',...
    '143319','143649','144058'};

expts(end+1).yymmdd = '20190813';
expts(end).hhmmss = {'113429','150813'};
expts(end).perf_hhmmss = {'114725','114930','151819','152724','153221','153437'};

% expts(end+1).yymmdd = '20190816';
% expts(end).hhmmss = {'110820'};
% expts(end).perf_hhmmss = {'112220','112435','112632','112926','113118','113218',...
%     '113506','114851','115108','115346'};

expts(end+1).yymmdd = '20190819';
expts(end).hhmmss = {'112813','141607'};
expts(end).perf_hhmmss = {'113941','114110','114318','114616','115111','115346',...
    '142226','142626','143335','143651','144113','144505','144816'};

expts(end+1).yymmdd = '20190830';
expts(end).hhmmss = {'112710','140606'};
expts(end).perf_hhmmss = {'114031','115045','115346','115737','141814','143231',...
    '143449','143810','144326','144640','145048','145735','150124'};

expts(end+1).yymmdd = '20190904';
expts(end).hhmmss = {'113405','140826'};
expts(end).perf_hhmmss = {'114718','115009','115356','115803','120024',...
    '142452','142649','142912','143715','143850','144117','144534','144902',...
    '145059','145426','145640','145824','150103'};


%% calculate and plot the convergence of KF in CLDA 

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
Expt_Ref=expts(end);
datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',Expt_Ref.yymmdd,...
    'GangulyServer','Center-Out',Expt_Ref.yymmdd,...
    Expt_Ref.hhmmss{1,length(Expt_Ref.hhmmss)},'BCI_CLDA');
files = dir(fullfile(datadir,'Data*.mat'));
load(fullfile(datadir,files(end).name));
C_Ref=TrialData.KalmanFilter{1,1}.C;

% calculate the angles compared to refs
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    % counting the days using CLDA sessions
    DayCLDASessions(i)=length(expt.hhmmss);
    
    for j=1:length(expt.hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\BRAVO1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.hhmmss{1,j},'BCI_CLDA');
        
        Sessions(end+1)={[yymmdd,'-',expt.hhmmss{1,j}]};
            
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

% Across sessions
%figure('position',[400 400 1800 300]);
figure('position',[445 326 1088 130]);
plot(Angle_x_AcrossSessions,'linewidth',1)
hold on
plot(Angle_y_AcrossSessions,'linewidth',1)
box off
set(gca,'FontSize',14)
xlim([0,max(cumsum(DayCLDASessions))])
ylim([0,90])
yticks(0:45:90)
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
hold on
vline(cumsum(DayCLDASessions))
title('KF weights convergence across all sessions')
%HighQualityFigsSVG('Fig_R1_1')

% Within sessions
figure('position',[400 400 1800 300]);
plot(Angle_x_WithinSessions,'linewidth',2)
hold on
plot(Angle_y_WithinSessions,'linewidth',2)
box off
set(gca,'FontSize',20)
xlim([0,max(cumsum(DayCLDASessions))])
ylim([0,90])
yticks(0:45:90)
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
hold on
vline(cumsum(DayCLDASessions))
title('KF weights variation within sessions')
%HighQualityFigsSVG('Fig_R1_2')

% Across trials
figure('position',[400 400 1800 300]);
plot(Angle_x_AcrossTrials,'linewidth',2)
hold on
plot(Angle_y_AcrossTrials,'linewidth',2)
box off
set(gca,'FontSize',20)
xticks(0:200:1000)
xlim([0,length(Angle_x_AcrossTrials)])
ylim([0,90])
yticks(0:45:90)
yt=get(gca,'ytick');
for k=1:numel(yt)
yt1{k}=sprintf('%d°',yt(k));
end
set(gca,'yticklabel',yt1);
hold on
vline(DayTrials)
title('KF weights convergence across all trials')
%HighQualityFigsSVG('Fig_R1_3')


%% Calculate the fitts rates, performance (success rates, time to target) corresponding to CLDA blocks

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    PerfsDays(i).TestDate=yymmdd;
    
    if i<12 %strcmp(TrialData.Params.Task,'Cetner-Out')
        FixedType='Center-Out';
        PerfsDays(i).FixedType=FixedType;
        
    else %strcmp(TrialData.Params.Task,'CursorControlGridTask')
        FixedType='CursorControlGridTask';
        PerfsDays(i).FixedType=FixedType;
    end 
       
    if i==1
       PerfsDays(i).TestDayCounter=0;
       Day_Ref=yymmdd;
    else 
      NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
          [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);
      PerfsDays(i).TestDayCounter=NumDays;
    end 
    
    % go througth the number of CLDA sessions
    for j=1:length(expt.hhmmss)
        fprintf('CLDA:')
        fprintf('\n%s-%s\n',yymmdd,expt.hhmmss{1,j})
        PerfsDays(i).CLDA(j).SessionNameCLDA=expt.hhmmss{1,j};
        
        if j<length(expt.hhmmss)
            
            % go througth the number of related fixed sessions
            FixedBlockCounter=0;
            for jj=1:length(expt.perf_hhmmss)
               
                if str2num(expt.perf_hhmmss{1,jj})>= str2num(expt.hhmmss{1,j}) && ...
                        str2num(expt.perf_hhmmss{1,jj})<str2num(expt.hhmmss{1,j+1})
                    
                    FixedBlockCounter=FixedBlockCounter+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_hhmmss{1,jj},FixedType);
                    
                    % save into structure:
                    PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).FittsRate=FittsRate;
                    PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).SuccessRate=SuccessRate;
                    PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).TravelTime=TravelTime;
                    PerfsDays(i).CLDA(j).FixedBlocks(FixedBlockCounter).Clicker=Clicker;
                                           
                end
                
            end
            
        elseif j==length(expt.hhmmss)
            
            FixedBlockCounter=0;
            for jj=1:length(expt.perf_hhmmss)
                
                if str2num(expt.perf_hhmmss{1,jj})>= str2num(expt.hhmmss{1,j})
                    
                    FixedBlockCounter=FixedBlockCounter+1;
                    %calculate the type of fixed block, fitts rate & performance & clicker
                    %on/off
                    [FittsRate, SuccessRate, TravelTime, Clicker]=FixedPerformanceFunc(yymmdd,expt.perf_hhmmss{1,jj},FixedType);
                    
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


% plotting the results
Perf_Days={};

for i=1:length(PerfsDays)
    
    Perf_Days(i)={num2str(PerfsDays(i).TestDayCounter)};
    counter=0;
    
    for j=1:length(PerfsDays(i).CLDA)
        
        for k=1:length(PerfsDays(i).CLDA(j).FixedBlocks)
            counter=counter+1;
            Perf_Plot(i).SuccessRate(counter)=PerfsDays(i).CLDA(j).FixedBlocks(k).SuccessRate;
            Perf_Plot(i).TravelTime(counter)=PerfsDays(i).CLDA(j).FixedBlocks(k).TravelTime;
            Perf_Plot(i).BitRate(counter)=PerfsDays(i).CLDA(j).FixedBlocks(k).FittsRate;
            Perf_Plot(i).Clicker(counter)={PerfsDays(i).CLDA(j).FixedBlocks(k).Clicker};
               
        end 
          
    end   
    
end 

% figure for marking the days for sessions_based plots
%figure('position',[400 400 1800 100]);
figure('position',[445 326 1088 70]);
box on
set(gca,'FontSize',12)
vline(cumsum(DayCLDASessions))
yticks('')
xlim([0,max(cumsum(DayCLDASessions))])
xticks(cumsum(DayCLDASessions))
xtickangle(45)
xticklabels(Perf_Days)
%HighQualityFigsSVG('Fig_R1_4')

% figure for bit rates
%figure('position',[400 400 1800 300]);
figure('position',[445 326 1088 110]);
SessionsNum=cumsum(DayCLDASessions);
for i=1:length(Perf_Plot)
    
    plot(SessionsNum(i),max(Perf_Plot(i).BitRate),'ok','MarkerFaceColor','k','MarkerSize',6)
    [m n]=max(Perf_Plot(i).BitRate);
    Index(i)=n;
    box off
    hold on  
end
yticks(0:0.4:0.8)
ylim([0 0.8])
xlim([0,max(cumsum(DayCLDASessions))])
set(gca,'FontSize',12)
hold on
vline(SessionsNum)
%HighQualityFigsSVG('Fig_R1_5')

% figure for success rate and travel time
%fig=figure('position',[400 400 1800 300]);
fig=figure('position',[445 326 1088 110]);
left_color = [1 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
for i=1:length(Perf_Plot)
    
    plot(SessionsNum(i),Perf_Plot(i).SuccessRate(Index(i)),'or','MarkerFaceColor','r','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,105])
xlim([0,max(cumsum(DayCLDASessions))])
set(gca,'FontSize',12)

yyaxis right
for i=1:length(Perf_Plot)
    
    plot(SessionsNum(i),Perf_Plot(i).TravelTime(Index(i)),'ob','MarkerFaceColor','b','MarkerSize',6)
    box off 
    hold on
end 
ylim([0,15])
set(gca,'FontSize',12)

hold on
%vline(SessionsNum)
%HighQualityFigsSVG('Fig_R1_6')








%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
