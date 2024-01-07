%% EE473 Final Project : LPC Vocoder
% Batuhan Tosun 2017401141

%% Overlap & Add Method 
% speech signal is divided into 30ms frames

[speech, fs] = audioread("my_input.wav");

speech = speech(:,1)';

% normalize speech
speech = speech/(max(abs(speech)));

% total number of samples for 30ms frame
N = floor(15*fs/1000); 

% total number of samples per window shift 15ms
R = floor(0.5*N);

% total number of frames created is equal to n_overlap 
n_overlap = floor((length(speech)-N)/R)+1;

% speech frames will be stored in X_m
X_m = zeros(n_overlap, N);

% Hanning window is used
w = hann(N, 'periodic');

for i=1:n_overlap
    
    X_m(i,:) = w'.*speech((1:N)+(i-1)*R);

end

%% Autocorrelation Method
% the filter coefficients for each frame are calculated using the autocorrelation method 

% total number of predictor coefficients (LP coefficients)
p = 32;

% matrix for holding the values of autocorrelation function of each frame
R_m = zeros(n_overlap,p+1);

% matrix for holding the short-time predictor coefficients for each frame
A = zeros(p,n_overlap);

% matrix for holding the error signal at the end of each frame
E = zeros(n_overlap,N);

for k = 1:n_overlap % for each frame
    
    s_frame = X_m(k,:);
    
    for j=1:p+1
        R_m(k,j) = s_frame(1:N-j+1)*s_frame(1+j-1:N)';
    end
    
    r1 = R_m(k,2:end);
    r2 = R_m(k,1:end-1);
    R_toep = toeplitz(r2);
    
    % alfa_k's are found
    A(:,k)=R_toep\r1';
    
    subConv = conv([0 X_m(k,:)], A(:,k)');
    E(k,:) = X_m(k,:)-subConv(1:N);
    
end


%% finding pitch periods (Center-clipping)

method1 = "center_clipping";
percentage_threshold=0.3;
[pitch_periods,pitch_frequencies] = pitch_estimation_func(X_m, fs, method1, percentage_threshold);

%% finding voiced & unvoiced (pitch-freq)

method = "pitch_frequency";
threshold_val = 200;
voiced_unvoiced = voiced_unvoiced_detector(X_m,method,threshold_val, pitch_frequencies);

%% gain calculation
gain = gain_calculator(R_m, A, pitch_periods,voiced_unvoiced);

%% Decoding section

Decoded_speech = zeros(n_overlap,N);

pulse_train = zeros(1,N);
offset = 0;

for i=1:n_overlap % for each speech frame
    
     % if the current speech frame is voiced we generate impulse train with a period equal to corresponding pitch period
    if(voiced_unvoiced(1,i)==1) 
        
        % creating impulse train
        src1 = zeros(N,1);
        step = round(pitch_periods(1,i)*fs);
        pts = (offset+1):step:N;
        
        if ~isempty(pts)
            offset = step + pts(end) - N;
            src1(pts) = 1; % impulse train, compensate power
        end
    
    % the resulting impulse train is multiplied by the corresponding gain
    % factor and then filtered out using the corresponding short-time
    % predictor coefficients 
        
        src1 = src1';
        pitch_periods(1,i)
        sum(src1)
        x_hat = filter(1,[1; -A(:,i)], src1*gain(1,i));
    
    % if the current speech frame is unvoiced we generate gaussian random
    % noise
    else
        src2 = randn(1,N);
        % the resulting gaussian is multiplied by the corresponding gain
        % factor and then filtered out using the corresponding short-time
        % predictor coefficients         
        x_hat = filter(1,[1; -A(:,i)],src2*gain(1,i));
    end
    % each decoded frame is stored in the rows of "Decoded_speech" matrix
    Decoded_speech(i,:)=x_hat;
end

%% decoding our speech

% here the decoded speech frames are concatenated in the reverse order of OLA and the decoded speech signal is created 
dec_speech = zeros(1,length(speech));
w = hann(N, 'periodic');

for j=1:n_overlap
    
    dec_speech = dec_speech + [zeros(1,(j-1)*R) w'.*Decoded_speech(j,:) zeros(1,length(speech)-(j-1)*R-N)];
    
end

% normalization
dec_speech = dec_speech/(max(abs(dec_speech)));

%% plotting

figure
subplot(2,1,1)
plot([0:length(speech)-1]/fs,speech)
title("Original Speech Signal")
xlabel("time(s)")
subplot(2,1,2)
plot([0:length(speech)-1]/fs,dec_speech)
title("LPC Decoded Speech Signal")
xlabel("time(s)")

figure
subplot(2,1,1)
plot([0:length(speech)-1]/fs,speech)
title("Original Speech Signal (1-2sec)")
xlabel("time(s)")
xlim([1,2])
subplot(2,1,2)
plot([0:length(speech)-1]/fs,dec_speech)
title("LPC Decoded Speech Signal (1-2sec)")
xlabel("time(s)")
xlim([1,2])

%%
disp("Original Sound Playing...")
sound(speech,fs)

pause(round(length(speech)/fs)+1);

disp("LPC-decoded Sound Playing...")
sound(-dec_speech,fs)
%%

audiowrite("my_output.wav",dec_speech,fs)