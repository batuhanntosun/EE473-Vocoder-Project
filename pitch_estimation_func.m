function [pitch_periods,pitch_frequencies] = pitch_estimation_func(X_m,fs, method, percentage_threshold)

% this function calculates the pitch period/frequency of each speech frame
% (chunk) using one of the two methods "autocorrelation" or "center clipping"
n_overlap = size(X_m,1);

pitch_periods = zeros(1,n_overlap);
pitch_frequencies = zeros(1,n_overlap);

for i = 1:n_overlap % for each speech frame
    
    if method == "autocorrelation"
        
        % we diractly take the autocorrelation function of the current
        % speech frame
        
        [R_xx, ~] = xcorr(X_m(i,:), X_m(i,:));
        % only the one side is enough
        R_xx_one_side = R_xx((length(R_xx)-1)/2+1:end);
        % then we find locations of the peaks (local maximum)
        [pks_max, locs_max] = findpeaks(R_xx_one_side);
        %the mean of the differences between local maximums is divided by
        %the sampling frequency and the result is our pitch period
        period_in_samples = mean(diff(locs_max));
        pitch_periods(1,i) = period_in_samples/fs;
        pitch_frequencies(1,i) = fs/period_in_samples;
       
    elseif method == "center_clipping"
        
        % first clipped version of the current frame is generated according
        % to the formula
        y_m = zeros(1,size(X_m,2));
        max_abs_val = max(abs(X_m(i,:)));
        clipping_threshold = percentage_threshold*max_abs_val;
        y_m = y_m+(X_m(i,:)-clipping_threshold).*(X_m(i,:)>=clipping_threshold);
        y_m = y_m + (X_m(i,:)+clipping_threshold).*(X_m(i,:)<=-clipping_threshold);
        
        % then the autocorrelation function of the center-clipped speech
        % frame is calculated
        
        [R_yy, ~] = xcorr(y_m, y_m);
        R_yy_one_side = R_yy((length(R_yy)-1)/2+1:end);
        %only the one side is enough
        % then we find locations of the peaks (local maximum)
        [pks_max, locs_max] = findpeaks(R_yy_one_side);
        %the mean of the differences between local maximums is divided by
        %the sampling frequency and the result is our pitch period        
        period_in_samples = locs_max(1);
        %period_in_samples = mean(diff(locs_max));
        pitch_periods(1,i) = period_in_samples/fs;
        pitch_frequencies(1,i) = fs/period_in_samples;    
        
    end

end

end

