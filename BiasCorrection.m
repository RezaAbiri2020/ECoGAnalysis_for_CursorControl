

cd('/home/ucsf/Projects/CursorControlGridTask/TaskCode/persistence')
load('kf_params.mat')


C=KF.C;
Q=KF.Q;
C = [C(129:256,:); C(769:end,:)];
Q1 = Q(129:256,129:256);
Q2 = Q(129:256,769:end);
Q3 = Q(769:end,129:256);
Q4 = Q(769:end,769:end);
Q = [Q1 Q2;Q3 Q4];

[C,Q,Qinv,tol,feat,angles]=bias_correct(C,Q,0,8,1,2.5);
[C,Q,Qinv,tol,feat,angles]=bias_correct(C,Q,180,8,1,2.5);


[KF.C,KF.Q,KF.Qinv,tol,feat,angles]=bias_correct(KF.C,KF.Q,0,8,1,2.5);
[KF.C,KF.Q,KF.Qinv,tol,feat,angles]=bias_correct(KF.C,KF.Q,180,8,1,2.5);

% to get rid of delta phase
% KF.C = KF.C(129:end,:);
% KF.Q = KF.Q(129:end,129:end);

FeatureMask(1:128) = 0;
FeatureMask = logical(FeatureMask);


%%%%% THIS IS OPTIONAL AND PROB BE UNSTABLE %%%%%%
% 
% % mask on the top 50% of features within each band 
% indices=[];
% for i=1:128:length(KF.Q) 
%     temp = KF.C(i:i+128-1,3) + 1i*KF.C(i:i+128-1,4);
%     [A I] = sort(abs(temp),'descend');
%     indices = [indices; I(1:64)];
% end
% 
% KF.C = KF.C(indices,:);
% KF.Q=KF.Q(indices,indices);
% % condition the matrix
% KF.Q = KF.Q + 1e-5*eye(size(KF.Q));
% KF.Qinv = inv(KF.Q);
% FeatureMask = FeatureMask(indices);

%%%%%%%%%%%% END %%%%%%%%

clearvars -except KF FeatureMask
save kf_params
