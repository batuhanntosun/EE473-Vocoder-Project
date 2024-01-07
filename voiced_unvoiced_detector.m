function [voiced_unvoiced] = voiced_unvoiced_detector(X_m, method, threshold_val, pitch_frequencies)

% this function classifies each speech frame as voiced or unvoiced
% (chunk) using one of the three methods: "zero-crossing",
% "pitch_frequency", "energy"

n_overlap =  size(X_m,1);
voiced_unvoiced = zeros(1,n_overlap);

for i = 1:n_overlap
    
    if method == "zero-crossing"
        
        signed_frame = sign(X_m(i,:));
        num_zero_cross = sum(abs(diff(signed_frame)));
        % if the total number of zero-crossing of the current speech frame is smaller than the
        % threshold value it is classified as "voiced", otherwise "unvoiced"
        voiced_unvoiced(1,i)= num_zero_cross<threshold_val;
    
    elseif method == "pitch_frequency"

        % if the estimated pitch frequency of the current speech frame is lower than the
        % threshold value it is classified as "voiced", otherwise "unvoiced"        
        voiced_unvoiced(1,i) = pitch_frequencies(1,i)<threshold_val;
           
    elseif method == "energy"
        % if the energy of the current speech frame is higher than the
        % threshold value it is classified as "voiced", otherwise "unvoiced"        
        energy_frame = sum(abs(X_m(i,:)).^2);
        voiced_unvoiced(1,i)= energy_frame>threshold_val;
    end
    
end    

end

