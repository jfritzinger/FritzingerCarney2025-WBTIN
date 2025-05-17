function [R2, avModel, stdModel, ratio, max_all] = modelWBTINSTRF(param_WB, data_STRF, data_WB)
% Predict neural responses using STRF model
%   [R2, AVG_MODEL, STD_MODEL, RATIO, MAX_ALL] = MODELWBTINSTRF(PARAM_WB, DATA_STRF, DATA_WB)
%   predicts WB-TIN responses by convolving stimuli with noise STRF.
%
%   Inputs:
%       param_WB - Struct containing stimulus parameters:
%           .physio : Physiological flag (0/1)
%           .Fs     : Sampling rate (Hz)
%           .stim   : Stimulus matrix
%           .mlist  : Metadata array with SNR/fpeak fields
%       data_STRF - STRF data struct:
%           .H2ex_strf : Excitatory STRF
%           .H2in_strf : Inhibitory STRF
%           .f         : Frequency vector (Hz)
%       data_WB   - Experimental data struct:
%           .rate : Measured response rates
%
%   Outputs:
%       R2        - R² values for each SNR level
%       avModel   - Predicted mean responses [fpeaks × SNRs]
%       stdModel  - Prediction standard deviations
%       ratio     - Experimental/model max response ratio
%       max_all   - Combined maximum response value
%

% Get STRF data 
STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
f = data_STRF.f;
f(f>16000) = [];
STRF_mat(length(f)+1:end, :) = [];

% Recreate stimulus
param_WB.physio = 1;
[param_WB] = generate_WBTIN(param_WB);

% Calculate the spectrogram of a stimulus
windowLength = 200; %50; % Length of the analysis window (in samples)
noverlap = 199; %49; % Overlap size between consecutive windows (in samples)
nfft = 960; % Length of the FFT

% Convolve STRF and the stimulus
num_stim = length(param_WB.list);
all_avg_rate = zeros(1, num_stim);
for istim = 1:num_stim

	% Plot the spectrogram of a stimulus
	[S,freq, ~] = spectrogram(param_WB.stim(istim, :), windowLength, noverlap, nfft, param_WB.Fs);
	freq(freq>16000) = [];
	S(length(freq)+1:end, :) = [];
	S = abs(S);

	% Convolve
	for ifreq = 1:length(freq)
		signal = S(ifreq,:);
		filter = STRF_mat(ifreq,:);
		conv1 = conv(signal, filter, 'same');
		convolution_result(ifreq,:) = conv1;
	end

	% Get PSTH & average rate
	psth = sum(convolution_result, 1);
	avg_rate = mean(psth);
	all_avg_rate(istim) = avg_rate;
end

[SNRs,~,si] = unique([param_WB.mlist.SNR].');
num_SNRs = length(SNRs);

[fpeaks,~,fi] = unique([param_WB.mlist.fpeak].');
num_fpeaks = length(fpeaks);

rate_size = [num_fpeaks,num_SNRs];
[avModel,stdModel,~,~] = accumstats({fi,si},all_avg_rate, rate_size);

% Model
noise_alone = mean(data_WB.rate(1,1:num_SNRs));
noise_alone_m = mean(avModel(1));
max_rate = max(max(data_WB.rate-noise_alone));
max_rate_m = max(max(avModel-noise_alone_m));
ratio =  max_rate/max_rate_m;
max_all = max([max_rate, max(((avModel-noise_alone_m).*ratio)+noise_alone)]);

% Measure prediction accuracy
for ir2 = 1:num_SNRs
	R_int = corrcoef(avModel(:, ir2),data_WB.rate(:,ir2));
	R2(ir2) = R_int(1, 2).^2;
end

end
