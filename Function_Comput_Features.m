
function [Output,RawData_trial]=Function_Comput_Features_V3(ECoG_data,Start)

%% this function will compute the seven featues from raw data collected from fixed blocks in center out task
Output={};

Burst_time1=0; % ms: start of the time for processing ERP
Burst_time2=2000; %ms: end of the time for processing ERP
Fs=1000;

% for just output
if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    RawData_trial={ECoG_data(((Start+Burst_time1)+1):((Start+Burst_time1)+1+Burst_time2)-1,:)};
    ECoG_data_1=cell2mat(RawData_trial);
else
    RawData_trial={[ECoG_data(((Start+Burst_time1)+1):end,:);NaN(Burst_time2-(size(ECoG_data,1)-(Start+Burst_time1)),size(ECoG_data,2))]};
    ECoG_data_1=ECoG_data(((Start+Burst_time1)+1):end,:);
end

% add white noise for edge artifact
m=4000;

%power=0;
%ECoG_data_2=[wgn(m,size(ECoG_data_1,2),power); ECoG_data_1;wgn(m,size(ECoG_data_1,2),power)]; 
Varation=std(ECoG_data_1);
ECoG_data_2=[repmat(Varation,m,1).*randn(m,size(ECoG_data_1,2)); ECoG_data_1;repmat(Varation,m,1).*randn(m,size(ECoG_data_1,2))]; 

%% Feature 1

[b,a]=butter(3,[0.5 4]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
Output_F1=exp(1i*angle(hilbert(ECoG_Filtered)));

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(1)={Output_F1((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(1)={[Output_F1((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end


%% Feature 2
[b,a]=butter(3,[0.5 4]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
Output_F2=abs(hilbert(ECoG_Filtered));

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(2)={Output_F2((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(2)={[Output_F2((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end



%% Feature 3
[b,a]=butter(3,[4 8]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
Output_F3=abs(hilbert(ECoG_Filtered));

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(3)={Output_F3((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(3)={[Output_F3((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end


%% Feature 4
[b,a]=butter(3,[8 13]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
Output_F4=abs(hilbert(ECoG_Filtered));

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(4)={Output_F4((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(4)={[Output_F4((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end


%% Feature 5
[b,a]=butter(3,[13 19]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
beta1=abs(hilbert(ECoG_Filtered));

[b,a]=butter(3,[19 30]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
beta2=abs(hilbert(ECoG_Filtered));

Output_F5=(beta1+beta2)/2;

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(5)={Output_F5((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(5)={[Output_F5((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end

%% Feature 6
[b,a]=butter(3,[30,36]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
lowgamma1=abs(hilbert(ECoG_Filtered));

[b,a]=butter(3,[36,42]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
lowgamma2=abs(hilbert(ECoG_Filtered));

[b,a]=butter(3,[42,50]/(Fs/2));
ECoG_Filtered=filtfilt(b,a,ECoG_data_2);
% hilbert
lowgamma3=abs(hilbert(ECoG_Filtered));

Output_F6=(lowgamma1+lowgamma2+lowgamma3)/3;

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(6)={Output_F6((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(6)={[Output_F6((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end


%% Feature 7
% pulling out Filtered data from different bands of HG
HGbands={
    [70,77]
    [77,85]
    [85,93]
    [93,102]
    [102,113]
    [113,124]
    [124,136]
    [136,150]};

for band=1:8
    [b,a]=butter(3,HGbands{band}/(Fs/2));
    FilteredAll(:,:,band)=filtfilt(b,a,ECoG_data_2);
    % Envelope & Hilbert
    [Upper_Envelope,Lower]=envelope(FilteredAll(:,:,band));
    % Hilbert & Filter
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HilbertAll(:,:,band)=filtfilt(b,a,Upper_Envelope)+repmat(mean(Upper_Envelope),size(Upper_Envelope,1),1);
    
end

Output_F7=mean(HilbertAll,3);

if (size(ECoG_data,1)-(Start+Burst_time1))>=Burst_time2
    Output(7)={Output_F7((m+1):(m+1+Burst_time2)-1,:)};
    
elseif (size(ECoG_data,1)-(Start+Burst_time1))<Burst_time2
    A=(size(ECoG_data,1)-(Start+Burst_time1));
    Output(7)={[Output_F7((m+1):(m+A),:);NaN(Burst_time2-A,size(ECoG_data,2))]};
    
end



end 