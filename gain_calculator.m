function [gain] = gain_calculator(R_m, A_m, pitch_periods,voiced_unvoiced)

% this function calculates the required filter gain for each speech frame
% according to the gain formula given in the lecture book
n_overlap =  size(R_m,1);
gain = zeros(1,n_overlap);

for i = 1:n_overlap % for each frame
    
    inner = R_m(i,1) - sum(A_m(:,i)'.*R_m(i,2:end));
    
    if voiced_unvoiced(1,i) == 1 % the formula changes wrt speech frame being voiced/unvoiced
        gain(1,i) = sqrt(pitch_periods(1,i)*inner);
    else
        gain(1,i) = sqrt(inner);
    end
    
end    

end

