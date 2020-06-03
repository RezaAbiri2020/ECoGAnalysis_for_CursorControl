function [Cnew,Qnew,Qinv,tol,feat,angles] =  bias_correct(C,Q,pref_angle,num_feat,bias_method,scale)
%
%function [Cnew,Qnew,Qinv,tol,feat,angles] =  bias_correct(C,Q,pref_angle,num_feat,bias_method,scale)
%
% Finds features across all freq bands and changes them to be biased
% towards a given direction, scaled to the max. feature in C matrix
%
% INPUT
% C             - C Matrix
% Q             - Q Matrix
% pref_angle    - Bias to be corrected. Angle measured clockwise, 0 deg
%                 being right, 90 deg being up, 180 deg being left, -180
%                 being left, -90 deg being down, -0 being right etc.
% num_feat      - number of features to be changed
% bias_method   - Whether to just change the gain of the feature (1) or
%                 scale it to the largest feature magnitide (2)
% scale         - If bias_method 1, then just the gain (e.g. 2). If
%                 bias_method 2 then between 0 and 1. i.e., 0.5 -> Scale to
%                 50% of largest magnitude in dataset
%
% OUTPUT
% Cnew - new C matrix
% Qnew - new Q matrix
% Qinv - new Qinv matrix
% tol - tolerance in deg (+/-) around the preferred angle. 
% feat - indices of features that were changed 
% angles - the angle differences between features and req. angle
%


Cnew = C;
Qnew = Q;

% convert preferred angle from degrees to radians
pref_angle = pref_angle * (pi/180);

% Get all weights
weights_all = C(1:end,3) + 1i*C(1:end,4);
weights_angle =  angle(weights_all);
weights_mag = abs(weights_all);

% pick the largest features biased towards requested angle
tol = 0;
n = 0;
while n<=num_feat
    angle_diff = angle(exp(1i*(weights_angle - pref_angle))); %angle differences
    feat =  find(abs(angle_diff)<tol);
    n=length(feat);
    tol = tol+(0.1*pi/180); % 0.1 degree increments
end
angles = angle_diff(feat)*(180/pi); % difference in terms of degrees

if bias_method == 1
    % multiply the gain to the feature magnitudes 
    weights_all_new=weights_all;
    for i=1:length(feat)        
        weights_all_new(feat(i)) = weights_all(feat(i)) * scale;
    end   
    
    % reduce the error in the Q matrix proportional to the gain 
    % lower bound is the lowest error in the model 
    [smallest_error smallest_error_feature] = min(diag(Q));
    Qd = diag(Q);
    target_error = (Qd(feat)/2);
    I = target_error > smallest_error;
    Qdf = Qd(feat);
    Qdf(I==1) = target_error(I==1);
    Qdf(I==0) = smallest_error;
    
    % update C matrix
    Cnew(:,3) = real(weights_all_new);
    Cnew(:,4) = imag(weights_all_new);
    
    % update Q matrix
    for i=1:length(feat)
        Qnew(feat(i),feat(i)) = Qdf(i);
    end    
    Qinv=inv(Qnew);
end


if bias_method == 2
    % find the largest feature weight
    [largest_weight largest_weight_feature] = max(abs(weights_all));
    target_weight = scale*largest_weight;
    
    
    % scale features in mag accordingly
    weights_all_new=weights_all;
    for i=1:length(feat)
        r = target_weight/weights_mag(feat(i));
        weights_all_new(feat(i)) = weights_all(feat(i)) * r;
    end
    
    % scale errors in Q matrix accordingly
    [smallest_error smallest_error_feature] = min(diag(Q));
    target_error = smallest_error + (1-scale)*smallest_error;
    
    % update C matrix
    Cnew(:,3) = real(weights_all_new);
    Cnew(:,4) = imag(weights_all_new);
    
    % update Q matrix
    for i=1:length(feat)
        Qnew(feat(i),feat(i)) = target_error;
    end
    Qinv=inv(Qnew);
end




