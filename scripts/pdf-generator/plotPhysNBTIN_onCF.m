function [params, fig, data] = plotPhysNBTIN_onCF(cluster, params, stims, CF, data)
% plotPhysNBTIN_onCF plots the on-CF average rate response to NB-TIN
%   This function finds the tone frequency nearest to the CF of the neuron
%   and plots the average rate response as a function of SNR.
%
% Inputs:
%    cluster - struct, output of post_process
%    params - params from post_process that have the type 'SPEC_slide' and
%    SPEC_slide_type 'NB_noise'
%	 stim - struct, output of post_process
%	 CF - characteristic frequency of the neuron
%	 data - cell array of structs for NB-TIN data
%
% Outputs:
%	 params - params from input, added plot_type to correctly plot PDFs.
%    fig - Average rate plot
%    data - cell of structs that contains processed average rate data.
% 		 data{ind}.rate_onCF - the on-CF average rates.
%		 data{ind}.p_onCF - p-value to see if there is a significant change
%		 between noise alone response and highest SNR response
%		 data{ind}.rate_std_onCF - standard deviation of average rate
%		 data{ind}.CF - CF of the neuron
%
% Other m-files required: accumstats, checkIfSPLShift
%
% Author: J. Fritzinger
% Created: 2023-07-24; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------

num_DSIDs = length(params);
binmode_strs = {'Ipsi','Contra','Bin'};
fig = figure;

for ind = 1:num_DSIDs

	% Parameters
	param = params{ind};
	ds = param.dsid;
	if isempty(stims)
		this_ds = param.stims.dsid == ds;
		cluster = param.cluster;
	else
		this_ds = stims(ds).dsid == ds;
	end
	binmode = binmode_strs{param.binmode + 1};

	% Analysis
	if param.version <= 5 % Not valid for plotting
		h = subplot(1, 1, 1);
		pos = get(h, 'position');
		dim = pos.*[1 1 0.5 0.5];
		str =  'This presentation was recorded before matching to WB-TIN';
		annotation('textbox',dim, 'String',str, 'FitBoxToText','on', 'verticalalignment', 'bottom')
		set(gca,'FontSize',7)
		params{1}.plot_type = 'NBTIN_onCF';
		params{1}.num_DSIDs = num_DSIDs;
		data{ind}.rate_onCF = [0 0 0];
		data{ind}.rate_std_onCF = [0 0 0];
		data{ind}.p_onCF = [];
		data{ind}.CF = [];
	elseif isnan(CF) % No CF recorded
		h = subplot(1, 1, 1);
		pos = get(h, 'position');
		dim = pos.*[1 1 0.5 0.5];
		str =  'No CF Recorded';
		annotation('textbox',dim, 'String',str, 'FitBoxToText','on', 'verticalalignment', 'bottom')
		set(gca,'FontSize',7)
		params{1}.plot_type = 'NBTIN_onCF';
		params{1}.num_DSIDs = num_DSIDs;
		box on
		title('On-CF NBTIN')
		data{ind}.rate_onCF = [0 0 0];
		data{ind}.rate_std_onCF = [0 0 0];
		data{ind}.p_onCF = [];
		data{ind}.CF = [];
	else
		% Check if this session has an spl shift
		spl_shift = checkIfSPLShift(param.animal, param.session);
		spl = param.spl + spl_shift;
		params{ind}.spl_real = spl;
		No = param.No + spl_shift;


		[fpeaks,~,fi] = unique([param.list.fpeak].');
		num_fpeaks = length(fpeaks);
		[SNR,~,si] = unique(double([param.list.SNR]).');
		SNR = flip(SNR);
		num_SNRs = length(SNR);

		if param.dur > 1
			dur = param.dur/1000; % stimulus duration in seconds.
		else
			dur = param.dur;
		end

		rate_size = [num_fpeaks,num_SNRs];
		spike_rates = cluster.num_spikes_delayed(this_ds)/...
			(dur - param.onsetWin/1000);

		if length(fi) == length(spike_rates)
			[rate,rate_std,~,~] = accumstats({fi,si},spike_rates, rate_size);

			SNRs(ind,:) = categorical(param.SNR);
			[~,closestIndex] = min(abs(fpeaks(:, 1)-CF));
			rate_onCF(ind, :) = rate(closestIndex,:);
			rate_std_onCF(ind,:) = rate_std(closestIndex,:);

			if sum(rate_onCF(ind,:))>0
				% T-test
				for j1 = 1:num_SNRs
					for j2 = closestIndex % bug in here
						y = spike_rates(fi == j2 & si == j1);
						rate_matrix(j1, 1:param.nrep) = y;
					end
				end
				[~,p] = ttest(rate_matrix(1,:), rate_matrix(3,:));
				if p<0.05
					label{ind} = sprintf('DS%d, No=%d, p = %.3f*', ds, No, p);
				else
					label{ind} = sprintf('DS%d, No=%d, p = %.3f', ds, No, p);
				end
			else
				label{ind} = 'Reanalyze, error';
				p = [];
			end

			data{ind}.rate_onCF = rate_onCF(ind, :);
			data{ind}.p_onCF = p;
			data{ind}.rate_std_onCF = rate_std_onCF(ind,:);
			data{ind}.CF = CF;
		else
			data{ind}.rate_onCF = [0 0 0];
			data{ind}.rate_std_onCF = [0 0 0];
			data{ind}.p_onCF = [];
			data{ind}.CF = [];
		end
	end
	clear leg
end

% Plot on-CF results
if param.version > 5 && ~isnan(CF)
	errorbar(SNRs', rate_onCF', (rate_std_onCF/sqrt(param.nrep))', 'Linewidth', 2);
	hLegend = legend(label, 'location', 'best', 'FontSize',6);
	hLegend.ItemTokenSize = [6,6];
	xlabel('SNR (dB SPL)')
	ylabel('Avg. Rate (sp/s')
	set(gca,'FontSize',7)
	grid on
	box on
	title(sprintf('On-CF NBTIN, CF=%dHz, Stim=%dHz', round(CF), round(fpeaks(closestIndex))));
	axis_set = ylim;
	axis_set(1) = 0;
	ylim(axis_set);
end

% Add type of plot
params{1}.plot_type = 'NBTIN_onCF';
params{1}.num_DSIDs = num_DSIDs;
end