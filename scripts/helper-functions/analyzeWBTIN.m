function data = analyzeWBTIN(params_WB, CF)
% Analyze wideband tone-in-noise (WBTIN) responses
%   data = analyzeWBTIN(params_WB, CF) processes neural responses to wideband
%   TIN stimuli with varying signal-to-noise ratios and frequency components.
%
%   Inputs:
%       params_WB - Cell array of parameter structures containing WB-TIN 
%					parameters and clustered data. 
%       CF - Characteristic frequency (Hz) for on-CF analysis. Enter [] to skip.
%
%   Outputs:
%       data - Cell array of structures containing analysis results:
%           Fields:
%           .rate           - Mean firing rate matrix (frequency × SNR)
%           .rate_std       - Standard deviation of firing rates
%           .rate_sm        - Smoothed rate vector (regularized fit)
%           .fpeaks         - Frequency components analyzed (Hz)
%           .V_p            - Predictable variance per SNR condition
%           .SNRs           - Signal-to-noise ratios analyzed (dB)
%           .rate_matrix    - Response matrix for variance calculation
%           .Asm            - Smoothed temporal response profile
%           .t_bins         - Time bin edges for temporal analysis (µs)
%           .rate_onCF      - CF-specific responses (if CF provided)
%           .rate_std_onCF  - CF response variability
%           .SNR_onCF       - SNR labels for CF responses
%           .ttest_sig      - p-value from CF response t-test
%
%   See also: accumstats, predictableVariance

% Process each parameter set
num_WB = length(params_WB);
data = cell(num_WB, 1);
for ind = 1:num_WB

	param = params_WB{ind};
	this_ds = param.stims.dsid==param.dsid;
	cluster = param.cluster;
	stim = param.stims;

	% Analysis
	[fpeaks,~,fi] = unique([param.list.fpeak].');
	num_fpeaks = length(fpeaks);
	[SNRs,~,si] = unique(double([param.list.SNR]).');
	num_SNRs = length(SNRs);

	if param.dur > 1
		dur = param.dur/1000; % stimulus duration in seconds.
	else
		dur = param.dur;
	end

	rate_size = [num_fpeaks,num_SNRs];
	fpeaks = repmat(fpeaks, 1, length(SNRs));
	spike_rates = cluster.num_spikes_delayed(this_ds)/...
		(dur - param.onsetWin/1000);

	if length(fi) == length(spike_rates)
		[rate,rate_std,~,rlb, rub] = accumstats({fi,si},spike_rates, rate_size);
		%rates_sm = smooth_rates(rate(:,2),rlb(:,2),rub(:,2), CF);

		rates_min = rlb(:,2);
		rates_max = rub(:,2);
		fun = @(x) sum(diff(diff(x)).^2) + 3*sum((x - rate(:,2)).^2);
		opt = optimset('fmincon');
		opt.Algorithm = 'active-set';
		opt.Display = 'off';
		rates_sm = fmincon(fun,rate(:,2),[],[],[],[],rates_min,rates_max,[],opt);

		data{ind}.rate_sm = rates_sm; 
		data{ind}.rate = rate;
		data{ind}.rate_std = rate_std;
		data{ind}.fpeaks = fpeaks;
	else
		data{ind}.rate = [];
		data{ind}.rate_std = [];
		data{ind}.fpeaks = [];
		data{ind}.V_p = [];
		data{ind}.SNRs = [];
		return
	end

	% On-CF Analysis
	if ~isempty(CF)
		
		% Find nearest frequency to CF
		[~, cfIdx] = min(abs(fpeaks(:,1) - CF));

		% Store baseline (noise alone) and CF responses
		data{ind}.rate_onCF = [mean(rate(1,:)), rate(cfIdx,:)];
		data{ind}.rate_std_onCF = [mean(rate_std(1,:)), rate_std(cfIdx,:)];
		data{ind}.SNR_onCF = categorical([-inf; param.SNR(:)], [-inf; param.SNR(:)]);

		% Compare lowest vs CF responses at highest SNR
		hiSNR = num_SNRs;
		baseRates = spike_rates(fi == 1 & si == hiSNR);
		cfRates = spike_rates(fi == cfIdx & si == hiSNR);

		% Ensure equal repetitions
		nrep = min(param.nrep, min(length(baseRates), length(cfRates)));
		[~, p] = ttest(baseRates(1:nrep), cfRates(1:nrep));

		data{ind}.ttest_sig = p;
	end

	% Calculate predictable variance per No
	for snr_ind = 1:num_SNRs
		rate_matrix = zeros([param.nrep, num_fpeaks]);
		x = reshape(spike_rates, [num_fpeaks*num_SNRs, param.nrep])';
		list_si = reshape(si, [num_fpeaks*num_SNRs, param.nrep])';
		list_fi = reshape(fi, [num_fpeaks*num_SNRs, param.nrep])';
		for irep = 1:param.nrep
			x_onerep = x(irep, :);
			si_onerep = list_si(irep, :);
			fi_onerep = list_fi(irep,:);
			for j1 = snr_ind
				for j2 = 1:num_fpeaks
					k = find(fi_onerep == j2 & si_onerep == j1);
					rate_matrix(irep, j2) = x_onerep(k);
				end
			end
		end
		[V_p(snr_ind), ~,~] = predictableVariance(rate_matrix, fpeaks);
	end
	data{ind}.V_p = V_p;
	data{ind}.SNRs = param.SNR;
	data{ind}.rate_matrix = rate_matrix;

	% Calculate port-fors
	noise_alone = rate(1,end);
	tinc = 12500; % size of the time bins, in µs
	tinc_s = tinc/1e6; % size of time bins, in seconds
	t_bins = 0:tinc:500000; % creates array of each time bin for the duration of stimulus
	num_bins = length(t_bins) - 1; % total number of bins
	portfors_bins = stim.times+t_bins; % Creates matrix of bin times
	portfors_bins = reshape(portfors_bins.',[],1); % Creates array of bin times

	% Determine number of spikes per 12.5 ms bin.
	num_portfors_spikes = histc(cluster.t_spike,portfors_bins);
	num_portfors_spikes(num_bins+1:num_bins+1:end) = [];
	num_portfors_spikes = reshape(num_portfors_spikes,num_bins,[]);

	ti = 1:num_bins;
	[ti2,fi2] = ndgrid(ti,fi);
	siz = [num_bins,num_fpeaks];

	SNR_ind = num_SNRs;
	k = si == SNR_ind;
	ti3 = ti2(:,k);
	fi3 = fi2(:,k);
	N3 = num_portfors_spikes(:,k);
	Tc = accumarray({ti3(:),fi3(:)},tinc_s,siz);
	%A = accumarray({ti3(:),fi3(:)},N3(:),siz)./Tc;
	A = accumarray({ti3(:),fi3(:)},N3(:),siz)./Tc - noise_alone;
	Asm = conv2(A,ones(3,3)/9,'same');
	data{ind}.Asm = Asm;
	data{ind}.t_bins = t_bins;

end
end
