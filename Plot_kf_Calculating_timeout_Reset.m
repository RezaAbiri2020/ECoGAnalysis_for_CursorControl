
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


%% Section 1-2: calculate Reset and timeout across sessions and days 

% finding the last trial as reference
Day_Ref='20190521';

Reset_Timeout=[];

for i=1:length(expts)
    expt = expts(i);
    yymmdd = expt.yymmdd;
        
    NumDays=daysdif([Day_Ref(1:4),'/',Day_Ref(5:6),'/',Day_Ref(7:8)],...
        [yymmdd(1:4),'/',yymmdd(5:6),'/',yymmdd(7:8)],1);    
    
    for j=1:length(expt.perf_hhmmss)
        fprintf('\n%s-%s\n',yymmdd,expt.perf_hhmmss{1,j})
        
        % go through datafiles in CLDA blocks
        datadir = fullfile('E:\Bravo1\CursorPlatform\Data\',yymmdd,...
            'GangulyServer','Center-Out',yymmdd,expt.perf_hhmmss{1,j},'BCI_Fixed');
        
        files = dir(fullfile(datadir,'Data*.mat'));
        load(fullfile(datadir,files(1).name));
        Reset_Timeout=[Reset_Timeout; ...
            NumDays,double(TrialData.Params.CenterReset),TrialData.Params.MaxReachTime];
        
        
    end
end


