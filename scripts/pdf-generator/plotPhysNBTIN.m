function [params, fig, data] = plotPhysNBTIN(cluster, params, stims, CF, type)
% plotPhysNBTIN plots the average rate response to sliding NB-TIN 
%   Average rate response for the NB-TIN stimulus from SPEC_slide.The
%   normal plot plots all average rates, the Subtract Noise plot subtracts
%   the noise-alone response from the other plots. 
%
% Inputs:
%    cluster - struct, output of post_process
%    params - params from post_process that have the type 'SPEC_slide' and
%    SPEC_slide_type 'NB_noise'
%	 stim - struct, output of post_process
%	 CF - characteristic frequency of the neuron
%	 type - 'Normal' or 'Subtract Noise'. 
%
% Outputs:
%	 params - params from input, added plot_type to correctly plot PDFs. 
%    fig - Average rate plot
%    data - cell of structs that contains processed average rate data. 
% 		 data{ind}.fpeaks - array of tone frequencies
%        data{ind}.rate - average rates
%        data{ind}.rate_std - standard deviation of average rates 
%		 data{ind}.SNRs - SNRs included in stimulus
%
% Other m-files required: accumstats, checkIfSPLShift
%
% Author: J. Fritzinger, adapted from cell_characterization
% Created: 2023-06-13; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------

num_DSIDs = length(params);
binmode_strs = {'Ipsi','Contra','Bin'};
binmode = binmode_strs{params{1}.binmode + 1};
fig = figure;
tiledlayout(ceil(num_DSIDs/3), 3, "Padding","compact", "TileSpacing","compact");
label_ind = 1;
data = cell(num_DSIDs, 1);

dsid_range_min = cellfun(@(p)p.fpeaks(1), params);
dsid_range_max = cellfun(@(p)p.fpeaks(end), params);
if length(unique(dsid_range_min)) == 1
    plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];
else
    [~,ID] = min(dsid_range_min);
    [~,ID_max] = max(dsid_range_max);
    plot_range = [params{ID}.fpeaks(1) params{ID_max}.fpeaks(end)];
end

for ind = 1:num_DSIDs
    param = params{ind};
    ds = param.dsid;
	if isempty(stims)
		this_ds = param.stims.dsid == ds;
		cluster = param.cluster;
	else
		this_ds = stims(ds).dsid == ds;
	end

    % Check if this session has an spl shift
    spl_shift = checkIfSPLShift(param.animal, param.session);
    spl = param.spl + spl_shift;
	params{ind}.spl_real = spl;

    % Analysis
    if ~isfield(param,'SNRNo') % Original NBTIN
        [fpeaks,~,fi] = unique([param.list.fpeak].');
        num_fpeaks = length(fpeaks);

        if param.dur > 1
            dur = param.dur/1000; % stimulus duration in seconds.
        else
            dur = param.dur;
        end

        rate_size = [num_fpeaks,1];
        spike_rates = cluster.num_spikes_delayed(this_ds)/...
            (dur - param.onsetWin/1000);

        if length(fi) == length(spike_rates)
            [rate,rate_std,~,~] = accumstats({fi},spike_rates, rate_size);

            % Plot
            hold on
            % Check this in older datasets! 
            errorbar(fpeaks,flip(rate, 2),flip(rate_std, 2)/(sqrt(param.nrep)));

            % Legend Labels
            if param.version == 1
                label(label_ind) = {sprintf('DSID %d, CF=%d', ...
                    param.dsid, param.fpeak_mid)};
            else
                label(label_ind) = {sprintf('DSID %d, CF=%d, SNR=%.2f', ...
                    param.dsid, param.fpeak_mid, param.SNR)};
            end
            label_ind = label_ind+1;

            % Save data 
            data{ind}.fpeaks = fpeaks;
            data{ind}.rate = rate;
            data{ind}.rate_std = rate_std;
        end
    elseif param.version >= 5 % Modified NBTIN to contain multiple SNRs in one dataset 
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

        nexttile
        if length(fi) == length(spike_rates)
            [rate,rate_std,~,~] = accumstats({fi,SNRi},spike_rates, rate_size);

            % Plot
            hold on
            %errorbar(fpeaks,flip(rate, 1),flip(rate_std,1)/(sqrt(param.nrep))); % Proper flip! 
			if strcmp(type, 'Normal')
				errorbar(fpeaks,rate,rate_std/(sqrt(param.nrep)));
			else
				errorbar(fpeaks,rate-rate(:,1),rate_std/(sqrt(param.nrep)));
			end

            % Legend Labels
            for ii = 1:num_SNRs
                leg_SNR{ii,:} = sprintf('%.1f SNR (%.0f)', SNRs(ii),param.SNRNo(ii));
            end
            plot_range = [param.fpeaks(1) param.fpeaks(end)];

            % Plot CF if known
            if ~isnan(CF)
                xline(CF, '--', 'Color', [0.5 0.5 0.5], 'LineWidth',1.5);
                leg_SNR{num_SNRs+1} = 'CF';
            else
                dim = [0.15 0.6 0.3 0.3];
                annotation('textbox',dim,'String','No CF stored','FitBoxToText','on');
            end
            hold off          

            % Label Figures
            set(gca, 'XScale', 'log');
            xlim(plot_range);
            xlabel('Frequency (kHz)')
			if strcmp(type, 'Normal')
				ylabel('Spike rate (sp/s)')
				axis_set = axis;
				axis_set(3) = 0;
				axis(axis_set);
			else
				ylabel('Spike rate - -inf rate (sp/s)')
			end
            set(gca,'FontSize',7)
            box on
            
            legend(leg_SNR, 'Location', 'best')
            if param.version == 5
                title('Tone not added!')
			else
				title(sprintf('CF=%d Hz, %.1f dB SPL, No=%d dB SPL, %s', ...
					param.fpeak_mid, spl, param.No, binmode));
			end

            % Save data
            data{ind}.fpeaks = fpeaks;
            data{ind}.rate = rate;
            data{ind}.rate_std = rate_std;
			data{ind}.SNRs = param.SNRNo;
        end
    end
end

if ~isfield(param,'SNRNo') % Original NBTIN
    if ~isnan(CF)
        xline(CF, '--', 'Color', [0.5 0.5 0.5]);
    else
        dim = [0.15 0.6 0.3 0.3];
        annotation('textbox',dim,'String','No CF stored','FitBoxToText','on');
    end
    label(label_ind) = {'CF'};
    hold off

    % Label Figures
    set(gca, 'XScale', 'log');
    xlim(plot_range);
    xlabel('Frequency (kHz)')
    ylabel('Spike rate (sp/s)')
    set(gca,'FontSize',7)
    grid on
    box on
    axis_set = axis;
    axis_set(3) = 0;
    axis(axis_set);
    hLegend = legend(label, 'Location', 'best');
	hLegend.ItemTokenSize = [6,6];

    title(sprintf('Third Octave Noise, %.1f dB SPL, %s, V.%d', ...
        spl, binmode,param.version));

    if param.version == 1 || param.version == 2
        dim = [0.15 0.6 0.3 0.3];
        annotation('textbox',dim,'String','Error in stimulus!!!','FitBoxToText','on');
    end
end

% Add type of plot
if strcmp(type, 'Normal')
	params{1}.plot_type = 'NBTIN';
else
	params{1}.plot_type = 'NBTIN_subnoise';
end
params{1}.num_DSIDs = num_DSIDs;

end