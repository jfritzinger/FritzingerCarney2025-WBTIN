function [params, fig, data] = plotPhysWBTIN(cluster, params, stims, CF, RM_rate, data, onset)
% Plots wideband tone-in-noise average rate response
% J. Fritzinger, updated 2/10/2022

num_DSIDs = length(params);
binmode_strs = {'Ipsi','Contra','Bin'};
fig = figure;
tiledlayout(ceil(num_DSIDs/3), 3, "Padding","compact", "TileSpacing","compact");
for ind = 1:num_DSIDs
	param = params{ind};
	ds = param.dsid;
	if isempty(stims)
		this_ds = param.stims.dsid == ds;
		cluster = param.cluster;
	else
		this_ds = stims(ds).dsid == ds;
	end
	binmode = binmode_strs{param.binmode + 1};

	% Check if this session has an spl shift
	spl_shift = checkIfSPLShift(param.animal, param.session);
	No = param.No + spl_shift;
	params{ind}.No_real = No;

	% Find BF nearest to WBTIN overall level
	lbw = param.fpeak_mid * 2^(-param.bandwidth/2);
	hbw = param.fpeak_mid * 2^(param.bandwidth/2);
	overall_spl = No + 10*log10(hbw-lbw);
	if ~isempty(RM_rate)
		[~,iBF] = min(abs(RM_rate.SPL-overall_spl));
		BF_spl = RM_rate.SPL(iBF);
		BF = RM_rate.BF(iBF);
	end

	% Analysis
	if param.version > 2
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

		if strcmp(onset, 'No')
			spike_rates = cluster.num_spikes_delayed(this_ds)/...
				(dur - param.onsetWin/1000);
		else
			spike_rates = cluster.num_spikes_peri(this_ds)/dur;
		end

		if length(fi) == length(spike_rates)
			[rate,rate_std,~,~] = accumstats({fi,si},spike_rates, rate_size);

			% Plot
			nexttile;
			hold on
			px = [fpeaks(2) fpeaks(end) fpeaks(end) fpeaks(2)];
			rate_pos = mean(rate(1,1:num_SNRs))+mean(rate_std(1,1:num_SNRs))/(sqrt(param.nrep));
			rate_neg = mean(rate(1,1:num_SNRs))-mean(rate_std(1,1:num_SNRs))/(sqrt(param.nrep));
			py = [rate_neg rate_neg rate_pos rate_pos];
			%patch(px, py, [0.9 0.9 0.9], 'EdgeColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3) % fix this
			patch(px, py, [0.9 0.9 0.9], 'EdgeColor', [0.8 0.8 0.8])
			line([fpeaks(2) fpeaks(end)], [1 1]* mean(rate(1,1:num_SNRs)),'Color','k');
			%errorbar(fpeaks,flip(rate, 1),flip(rate_std, 1)/(sqrt(param.nrep)));
			errorbar(fpeaks,rate,rate_std/(sqrt(param.nrep)));

			if ~isnan(CF)
				xline(CF, '--', 'Color', [0.5 0.5 0.5], 'LineWidth',1.5);
			else
				dim = [0.15 0.6 0.3 0.3];
				annotation('textbox',dim,'String','No CF value stored',...
					'FitBoxToText','on', 'FontSize',7);
			end

			if ~isempty(RM_rate)
				xline(BF, '--', 'Color', [1 0 0], 'LineWidth',1.5);
			end

			% Find Peak
			[max_rate, iMax] = max(rate, [],'all');
			max_rate_freq = fpeaks(iMax);
			scatter(max_rate_freq, max_rate, 'filled');

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

			% 			% Calculate overall predictable variance
			% 			rate_matrix = zeros(param.nrep, num_fpeaks,num_SNRs);
			% 			x = reshape(spike_rates, [num_fpeaks*num_SNRs, param.nrep])';
			% 			list_si = reshape(si, [num_fpeaks*num_SNRs, param.nrep])';
			% 			list_fi = reshape(fi, [num_fpeaks*num_SNRs, param.nrep])';
			% 			for irep = 1:param.nrep
			% 				x_onerep = x(irep, :);
			% 				si_onerep = list_si(irep, :);
			% 				fi_onerep = list_fi(irep,:);
			% 				for j1 = 1:num_SNRs
			% 					for j2 = 1:num_fpeaks
			% 						k = find(fi_onerep == j2 & si_onerep == j1);
			% 						rate_matrix(irep, j2, j1) = x_onerep(k);
			% 					end
			% 				end
			% 			end
			% 			a = reshape(rate_matrix,param.nrep, num_fpeaks*num_SNRs,1); % reshape
			% 			f = reshape(fpeaks, num_SNRs*num_fpeaks, 1);
			% 			[V_o, V_o2,~] = predictableVariance(a, f);


			hold off
			if param.version < 4
				lbw = param.fpeak_mid * 2^(-param.bandwidth/2);
				hbw = param.fpeak_mid * 2^(param.bandwidth/2);
				overall_level = No + 10*log10(hbw-lbw);
				new_SNRs = overall_level + SNRs - param.No;
				format short g
				for ii = 1:num_SNRs
					leg(ii,:) = [num2str(round(new_SNRs(ii),1)) 'dB SNR'];
				end
				hLegend = legend(char('Noise Alone', leg), 'location', 'best');
				hLegend.ItemTokenSize = [6,6];
			else
				for ii = 1:num_SNRs
					leg{ii} = sprintf('%2g dB SNR, Vp=%.2g', SNRs(ii), V_p(ii));
				end
				if ~isnan(CF) && ~isempty(RM_rate)
					BF_leg =  sprintf('BF at %ddB SPL', round(BF_spl));
					CF_leg = sprintf('CF');
					hLegend = legend(char('Noise Alone Error', 'Noise Alone',...
						leg{:},  CF_leg, BF_leg), 'location', 'best', 'FontSize',5);
					hLegend.ItemTokenSize = [6,6];
				elseif ~isnan(CF)
					CF_leg = sprintf('CF');
					hLegend = legend(char('Noise Alone Error', 'Noise Alone',...
						leg{:}, CF_leg), 'location', 'best', 'FontSize',5);
					hLegend.ItemTokenSize = [6,6];
				elseif ~isempty(RM_rate)
					BF_leg =  sprintf('BF at %d dB SPL', round(BF_spl));
					hLegend = legend(char('Noise Alone Error', 'Noise Alone',...
						leg{:}, BF_leg), 'location', 'best', 'FontSize',5);
					hLegend.ItemTokenSize = [6,6];
				else
					hLegend = legend(char('Noise Alone Error', 'Noise Alone', ...
						leg{:}), 'location', 'best', 'FontSize',5);
					hLegend.ItemTokenSize = [6,6];
				end
			end

			% Label Figures
			set(gca, 'XScale', 'log');
			if param.version >= 5
				plot_range = [param.fpeak_mid * 2^(-param.range/2), ...
					param.fpeak_mid * 2^(param.range/2)];
				if param.range > 3
					mid = param.fpeak_mid;
					x_ticks = [mid*2^-2, mid*2^-1 mid mid*2^1 mid*2^2];
					xticks(x_ticks)
				end
			else
				plot_range = [param.fpeak_mid/2 param.fpeak_mid*2];
			end
			xlim(plot_range);
			xlabel('Frequency (kHz)')
			ylabel('Spike rate (sp/s)')
			set(gca,'FontSize',7)
			box on
			axis_set = axis;
			axis_set(3) = 0;
			axis(axis_set);

			if strcmp(onset, 'No')
				title(sprintf('DS%d, CF=%d, No=%d, BW=%d, %s, SPL=%d', ...
					ds, param.fpeak_mid, No, param.bandwidth, ...
					binmode, round(overall_spl)));
			else
				title(sprintf('DS%d, No=%d, BW=%d, %s, With Onset', ...
					ds, No, param.bandwidth, ...
					binmode));
			end

			% Create struct that contains processed data
			data{ind}.fpeaks = fpeaks;
			data{ind}.rate = rate;
			data{ind}.rate_std = rate_std;
			data{ind}.spike_rates = spike_rates;
			data{ind}.V_p = V_p;
			data{ind}.SNRs = SNRs;
		else
			data{ind}.fpeaks = [];
			data{ind}.rate = [];
			data{ind}.rate_std = [];
			data{ind}.spike_rates = [];
			data{ind}.V_p = [];
		end
	else
		data{ind}.fpeaks = [];
		data{ind}.rate = [];
		data{ind}.rate_std = [];
		data{ind}.spike_rates = [];
		data{ind}.V_p = [];
		h = subplot(round(num_DSIDs/2), 2, ind);
		pos = get(h, 'position');
		dim = pos.*[1 1 0.5 0.5];
		str =  ['This presentation was version ' num2str(param.version)];
		annotation('textbox',dim, 'String',str, 'FitBoxToText','on', 'verticalalignment', 'bottom')
		set(gca,'FontSize',7)
	end
	clear leg
end

% Add type of plot
if strcmp(onset, 'No')
	params{1}.plot_type = 'WBTIN';
else
	params{1}.plot_type = 'WBTIN_Onset';
end
params{1}.num_DSIDs = num_DSIDs;

end
