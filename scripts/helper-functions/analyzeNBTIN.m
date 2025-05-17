function data = analyzeNBTIN(params, CF)
% Analyze narrowband TIN (tone-in-noise) responses
%   data = analyzeNBTIN(params, CF) processes neural responses to narrowband
%   stimuli with varying center frequencies and signal-to-noise ratios and
%   gives average rate response output. 
%
%   Inputs:
%       params - Cell array of parameter structures
%       CF - Characteristic frequency (Hz) for on-CF analysis. Enter [] to skip.
%
%   Outputs:
%       data - Cell array of structures containing analysis results for each parameter set
%           Fields:
%           .fpeaks          - Center frequencies analyzed (Hz)
%           .rate            - Mean firing rate matrix (frequency × SNR)
%           .rate_std        - Standard deviation of firing rates
%           .SNRs            - Signal-to-noise ratios analyzed (dB)
%           .rate_onCF       - Firing rates at CF (if CF provided)
%           .rate_std_onCF   - Standard deviation at CF
%           .SNR             - Categorical SNR labels
%           .spl             - Sound pressure level at CF (dB SPL)
%           .ttest_sig       - p-value from paired t-test (CF responses)
%
%   See also: accumstats

num_NB = length(params);
data = cell(num_NB, 1);
for ind = 1:num_NB
	param = params{ind};
	this_ds = param.stims.dsid==param.dsid;
	cluster = param.cluster;

	[fpeaks,~,fi] = unique([param.list.fpeak].');
	num_fpeaks = length(fpeaks);
	[SNRs,~,SNRi] = unique([param.list.SNR].');
	num_SNRs = length(SNRs);

	if param.dur > 1
		dur = param.dur/1000; % stimulus duration in seconds.
	else
		dur = param.dur;
	end

	rate_size = [num_fpeaks,num_SNRs];
	spike_rates = cluster.num_spikes_delayed(this_ds)/...
		(dur - param.onsetWin/1000);
	fpeaks = repmat(fpeaks, 1, length(SNRs));

	if length(fi) == length(spike_rates)
		[rate,rate_std,~,~] = accumstats({fi,SNRi},spike_rates, rate_size);
	end

	% Save data
	data{ind}.fpeaks = fpeaks;
	data{ind}.rate = rate;
	data{ind}.rate_std = rate_std;
	data{ind}.SNRs = param.SNRNo;

	% On-CF Analysis
	if ~isempty(CF)
		
		% Find frequency in fpeaks nearest to CF
		[~,closestIndex] = min(abs(fpeaks(:, 1)-CF));

		% Find rate nearest to CF
		SNRs(ind,1:num_SNRs) = categorical(param.SNRNo);
		rate_onCF(ind, :) = rate(closestIndex,:);
		rate_std_onCF(ind,:) = rate_std(closestIndex,:);

		% Find sound level at CF
		lbw = fpeaks(closestIndex) * 2^(-1/6);
		hbw = fpeaks(closestIndex) * 2^(1/6);
		spl = param.No + 10*log10(hbw-lbw);

		% T-test
		for j1 = 1:num_SNRs
			for j2 = closestIndex 
				y = spike_rates(fi == j2 & SNRi == j1);
				rate_matrix(j1, 1:param.nrep) = y;
			end
		end
		[~,p] = ttest(rate_matrix(1,:), rate_matrix(3,:));

		% Save data
		data{ind}.SNRs = param.SNRNo;
		data{ind}.rate_onCF = rate_onCF;
		data{ind}.rate_std_onCF = rate_std_onCF;
		data{ind}.SNR = SNRs;
		data{ind}.spl = spl;
		data{ind}.ttest_sig = p;
	end
end
end
