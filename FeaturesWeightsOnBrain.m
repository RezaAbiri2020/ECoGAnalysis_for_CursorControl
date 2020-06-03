

%% Plotting weights on the brain
% if all features zscored together or normalized together; then we can compare their weights on brain

clear all
close all
NCh=128;
%% loading fixed parameters for fixed block
% for first run: 5hz/5hz for fixed RML DC 
datadir = uigetdir();
datafiles = dir(fullfile(datadir,'Data*.mat'));

% load first trial
load(fullfile(datadir,datafiles(end).name))

%number of features
NFeatures=TrialData.Params.NumFeatures;
%Check if zscore happened for all features
zfeatures=TrialData.Params.ZscoreFeaturesFlag;

Coeff_Features=TrialData.KalmanFilter{1,end}.C;

% third column for vx and forth column for vy
%plot(Coeff_Features(:,3))
%plot(Coeff_Features(:,4))
 

% for vx
Vx_Coeffs=reshape(Coeff_Features(:,3),NCh,NFeatures);
% for vy
Vy_Coeffs=reshape(Coeff_Features(:,4),NCh,NFeatures);

% saving 
Vx_Coeffs_3=Vx_Coeffs;
Vy_Coeffs_3=Vy_Coeffs;
save('Apr17MorningCoeffs.mat','Vx_Coeffs_3','Vy_Coeffs_3')

%% Loading all three session data and analysis
load('Apr26MorningCoeffs')
load('Apr17MorningCoeffs')
load('May14AfternoonCoeffs')

Vx_Coeffs_All={Vx_Coeffs_3,Vx_Coeffs_1,Vx_Coeffs_2,(Vx_Coeffs_3+Vx_Coeffs_2+Vx_Coeffs_1)/3,Vx_Coeffs_3-Vx_Coeffs_1...
    ,Vx_Coeffs_3-Vx_Coeffs_2,Vx_Coeffs_1-Vx_Coeffs_2};

Vy_Coeffs_All={Vy_Coeffs_3,Vy_Coeffs_1,Vy_Coeffs_2,(Vy_Coeffs_3+Vy_Coeffs_2+Vy_Coeffs_1)/3,Vy_Coeffs_3-Vy_Coeffs_1...
    ,Vy_Coeffs_3-Vy_Coeffs_2,Vy_Coeffs_1-Vy_Coeffs_2};

%Feature 1: Phase of delta
%Feature 2: delta power
%Feature 3: theta power
%Feature 4: alpha power
%Feature 5: beta  power
%Feature 6: gamma1 power
%Feature 7: gamma2 power

NFeatures=7;
FeatureName={'phase of delta','delta power','theta power','alpha power','beta power','gamma1 power','gamma2 power'};
SessionTitles={'Apr17','Apr26','May14','Averaged Weigths','Apr17-Apr26','Apr17-May14','Apr26-May14'};

for i=1:NFeatures
    
    % for Vx
    figure(i),
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{i},' in x direction: Weights of the feature during fixed blocks'])
    
    % for Vy
    figure(i+NFeatures),
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle([FeatureName{i},' in y direction: Weights of the feature during fixed blocks'])
    
    
    for NCoeff=1:7
        
        Vx_Coeffs=Vx_Coeffs_All{NCoeff};
        Vy_Coeffs=Vy_Coeffs_All{NCoeff};
        % converting for correct mapping between pedestal and grid
        ch_layout = [
            91	84	67	90	70	79	88	69	92	83	65	89	87	86	94	82
            66	93	78	95	76	75	85	73	68	80	74	72	96	71	77	81
            60	37	42	50	56	54	49	40	43	35	45	63	47	46	58	55
            53	57	33	48	39	51	41	34	64	52	62	38	36	44	61	59
            8	26	29	28	9	5	13	20	11	23	16	22	27	4	3	31
            7	21	15	24	25	1	2	32	14	12	30	19	18	17	6	10
            110	125	111	115	103	117	100	123	113	119	118	98	101	105	116	99
            107	112	97	128	121	124	108	109	127	126	106	122	114	120	104	102];
        
        ch_layout=ch_layout';
        ch_layout=ch_layout(:);
        Vx_Coeffs_Brain=zeros(NCh,NFeatures);
        Vy_Coeffs_Brain=zeros(NCh,NFeatures);
        
        for j=1:NCh
            Vx_Coeffs_Brain(j,:)=Vx_Coeffs(ch_layout(j),:);
            Vy_Coeffs_Brain(j,:)=Vy_Coeffs(ch_layout(j),:);
        end
        
        % I have Coeffs_Brain for this structure for each feature:
        % 1 2 3 ..
        % 17 18....
        % ........128
        
        % I Should do numbering based on following structure for feeding into brain plots
        % codes
        % 113.....128
        %     .......
        % 1 2 3  ...16
        Vx_Coeffs_Brain_Mo=zeros(NCh,NFeatures);
        Vy_Coeffs_Brain_Mo=zeros(NCh,NFeatures);
        
        for ii=1:NFeatures
            A=reshape(Vx_Coeffs_Brain(:,ii),16,8);
            A=A';
            A=flip(A);
            A=A';
            A=A(:);
            Vx_Coeffs_Brain_Mo(:,ii)=A;
            
            B=reshape(Vy_Coeffs_Brain(:,ii),16,8);
            B=B';
            B=flip(B);
            B=B';
            B=B(:);
            Vy_Coeffs_Brain_Mo(:,ii)=B;
            
        end
        
        %% plot the coefficient of C on the brain
        
        load('BRAVO1_lh_pial')
        load('elecs_all')
        
        % for vx
        figure(i),
        subplot(2,4,NCoeff)
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Vx_Coeffs_Brain_Mo(:,i),'lh');
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        colorbar
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (SessionTitles(NCoeff))
       
        % for vy
        figure(i+NFeatures),
        subplot(2,4,NCoeff)
        ctmr_gauss_plot(cortex,elecmatrix(1:NCh,:),Vy_Coeffs_Brain_Mo(:,i),'lh');
        el_add(elecmatrix(1:NCh,:),'msize',1.7);
        colorbar
        %el_add(elecmatrix(1:128,:),'msize',1.7,'color', 'b', 'numbers', ch);
        title (SessionTitles(NCoeff))
        
        
        
        
    end
    
    figure(i),
    HighQualityFigs([FeatureName{i},'_WeightsVx'])
    figure(i+NFeatures),
    HighQualityFigs([FeatureName{i},'_WeightsVy'])
end
