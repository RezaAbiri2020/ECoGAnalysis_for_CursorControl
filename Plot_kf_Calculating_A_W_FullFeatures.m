
%% Section 1:calculating W and A using collected data to double check the vales: using fixed blocks

% should use fixed blocks
clear all 
close all
clc

%% Section 1-1: loading the data

%experiment info
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


%% Section 1-2: calculate W and A for trials across sessions and days 

axx_estim=[];
ayy_estim=[];
wxx_estim=[];
wyy_estim=[];
a_real=[];
w_real=[];
                
for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.perf_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.perf_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.perf_hhmmss{1,j},'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            
              if length(TrialData.Events)==2 % it means the cursor start from center 
                idx=find(TrialData.Time>=TrialData.Events(2).Time);
                Kinematics= TrialData.CursorState(:,idx);
                Gain=TrialData.Params.Gain;
                
                X1=[Gain*Kinematics(1:2,1:end-1);Kinematics(3:4,1:end-1);Kinematics(5,1:end-1)];
                X2=[Kinematics(1:2,2:end);Gain*Kinematics(3:4,2:end);Kinematics(5,2:end)];
                A=X2*X1'*pinv(X1*X1');
                
                % constrain A; first method
                if 0
                A_Ideal=zeros(5,5);
                A_Ideal(1,[1,2])=A(1,[1,2]);
                A_Ideal(2,[2,4])=A(2,[2,4]);
                A_Ideal(3,[3,4])=A(3,[3,4]);
                A_Ideal(4,[3,4])=A(4,[3,4]);
                A_Ideal(5,5)=A(5,5);
                A=A_Ideal;
                end 
                % constrain A; second method
                %A=[A(3,3),A(3,4);A(4,3),A(4,4)];
                %Gain=TrialData.Params.Gain;
                %A(1,3)=1/Gain*A(1,3);
                %A(2,4)=1/Gain*A(2,4);
                
                Value=(X2-A*X1)';
                W=1/(length(Kinematics)-1)*(X2-A*X1)*Value;
                
                axx_estim(end+1)=A(1,1);
                ayy_estim(end+1)=A(2,2);
                wxx_estim(end+1)=W(1,1);
                wyy_estim(end+1)=W(2,2);
                
                a_real(end+1)=TrialData.Params.KF.A(3,3);
                w_real(end+1)=TrialData.Params.KF.W(3,3);
                
              end       
        end   
    end
end

% G=?;
% t = 1/Params.UpdateRate;
% a=?
%  A = [...
%         1	0	G*t	0	0;
%         0	1	0	G*t	0;
%         0	0	a	0	0;
%         0	0	0	a	0;
%         0	0	0	0	1];
% w=?
% W = [...
%         0	0	0	0	0;
%         0	0	0	0	0;
%         0	0	w	0	0;
%         0	0	0	w	0;
%         0	0	0	0	0];

%% Section 1-3: plotting 
% for a 
figure('position',[445 326 1300 400]);
plot(a_real,'linewidth',1)
hold on
plot(axx_estim,'linewidth',1)
hold on 
plot(ayy_estim,'linewidth',1)
%ylim([0 1])
legend('a Real','axx Estim','ayy Estim')
xlabel('Trials')
title('a values in A matrix')
%HighQualityFigs('Matrix_A_variation_V1')

% for w 
figure('position',[445 326 1300 400]);
plot(w_real,'linewidth',1)
hold on
plot(wxx_estim,'linewidth',1)
hold on 
plot(wyy_estim,'linewidth',1)
%ylim([0 202])
legend('w Real','wxx Estim','wyy Estim')
xlabel('Trials')
title('w values in W matrix')
%HighQualityFigs('Matrix_W_variation_V2')

%% Section 2:calculating W and A using collected data to double check the vales: using only imagined blocks

% should use fixed blocks
clear all 
close all
clc

%% Section 2-1: loading the data

%experiment info
% imagined in Center-Out 
expts = [];

expts(end+1).yymmdd = '20190501';
expts(end).Imagined_hhmmss = {'133745'};

expts(end+1).yymmdd = '20190507';
expts(end).Imagined_hhmmss = {'105331'};

expts(end+1).yymmdd = '20190510';
expts(end).Imagined_hhmmss = {'104601','132311'};

expts(end+1).yymmdd = '20190514';
expts(end).Imagined_hhmmss = {'111822'};

expts(end+1).yymmdd = '20190515';
expts(end).Imagined_hhmmss = {'105804','132722'};



%% Section 2-2: calculate W and A for trials across sessions and days 

AllKinematics=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
    
    for j=1:length(expt.Imagined_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.Imagined_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.Imagined_hhmmss{1,j},'Imagined');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        for k=1:length(files)
            load(fullfile(datadir,files(k).name));
            idx=find(TrialData.Time>=TrialData.Events(2).Time);
            Kinematics= TrialData.CursorState(:,idx);         
            AllKinematics=[AllKinematics,Kinematics];
        end
    end
end

X1=AllKinematics(:,1:end-1);
X2=AllKinematics(:,2:end);
A=X2*X1'*pinv(X1*X1');

% constrain A; first method
if 0
    A_Ideal=zeros(5,5);
    A_Ideal(1,[1,2])=A(1,[1,2]);
    A_Ideal(2,[2,4])=A(2,[2,4]);
    A_Ideal(3,[3,4])=A(3,[3,4]);
    A_Ideal(4,[3,4])=A(4,[3,4]);
    A_Ideal(5,5)=A(5,5);
    A=A_Ideal;
end
% constrain A; second method

Value=(X2-A*X1)';
W=1/(length(AllKinematics)-1)*(X2-A*X1)*Value;

axx_estim=A(3,3);
ayy_estim=A(4,4);
wxx_estim=W(3,3);
wyy_estim=W(4,4);

% W is around 120 A is around 0.9

% G=?;
% t = 1/Params.UpdateRate;
% a=?
%  A = [...
%         1	0	G*t	0	0;
%         0	1	0	G*t	0;
%         0	0	a	0	0;
%         0	0	0	a	0;
%         0	0	0	0	1];
% w=?
% W = [...
%         0	0	0	0	0;
%         0	0	0	0	0;
%         0	0	w	0	0;
%         0	0	0	w	0;
%         0	0	0	0	0];

%% Section 2-3: plotting 



