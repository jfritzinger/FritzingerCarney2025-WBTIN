function [BMF, WMF, MTF_shape,fig, data] = plotPhysMTF_TTest(cluster, params, stims)
% plotPhysMTF_Ttest plots the average rate response to modulated noise.
%   Average rate response for amplitude modulated broadband gaussian noise.
%   The amplitude modulation varies from unmodulated to approximately 500
%   Hz
%
% Inputs:
%    cluster - struct, output of post_process
%    params - params from post_process that have the type 'typMTFN'
%	 stim - struct, output of post_process
%
% Outputs:
%	 BMF - Best modulation frequency
%	 WMF - Worst modulation frequency
%	 MTF_shape - MTF shape, BE, BS, H (BE-BS), H (BS-BE) or F
%    fig - Average rate plot
%    data - struct that contains processed average rate data. 
% 		 data.rate - average rate response
%        data.fms - modulation frequencies
%        data.rate_sm - smoothed average rate response
%        data.rate_std - standard deviation of rate response
%		 data.t_spike - from cluster, spike times 
%		 data.t_spike_rel - from cluster, spike times relative to dataset
%		 data.stim_times - from stim, times of stimulus
%		 data.abs_stim_num - from cluster
%		 data.at_100 - MTF shape at 100 Hz modulation
%		 data.MTF_shape - MTF shape, BE, BS, H (BE-BS), H (BS-BE) or F
%		 data.BMF - Best modulation frequency
%		 data.WMF - Worst modulation frequency
%
% Other m-files required: accumstats, checkIfSPLShift, MTFclassification,
%	smooth_rates 
%
% Author: J. Fritzinger, adapted from cell_characterization
% Created: 2023-01-30; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------

% Process data 
if iscell(params)
	param = params{1}; % If multiple, plots only the first instance
else
	param = params;
end
ds = param.dsid; 
if length(stims)==1
	this_ds = stims.dsid == ds;
	stim = stims;
else
    dsi = 1;
	this_ds = stims(ds).dsid == ds(dsi);
	stim = stims(ds);
end

fig = figure;

% Check if this session has an spl shift 
spl_shift = checkIfSPLShift(param.animal, param.session);
No = param.No +spl_shift;

if param.dur > 10
    dur = param.dur/1000; % stimulus duration in seconds.
else
    dur = param.dur;
end

all_mod_depths = double([param.list.mdepth]).';
all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
[~,~,mdi] = unique(all_mod_depths);
[fms,~,fmi] = unique(double([param.list.fm]).');
if fms(1) == 0
    fms(1) = 1.2;
end

num_mod_freqs = length(fms);
num_depths = length(param.all_mdepths);
map_size = [num_mod_freqs, num_depths];

spike_rates = cluster.num_spikes_delayed(this_ds)/...
    (dur - param.onsetWin/1000);

if length(fmi)==length(spike_rates)
    [rate,rate_std,rlb,rub] = accumstats({fmi,mdi},spike_rates, map_size);
    rate_sm = zeros(size(rate));
    for j = 1:num_depths
        rate_sm(:,j) = smooth_rates(rate(:,j),rlb(:,j),rub(:,j), []);
    end

    for j = 1:num_mod_freqs
        for k = 1:num_depths
            indices = j==fmi & k == mdi;
            raw_rates(j, k, :) = spike_rates(indices);
        end
    end

    % Plot
    ax = axes;
    binmode_strs = {'Ipsi','Contra','Bin'};
    if param.binmode == 8
        binmode = 3;
    else
        binmode = binmode_strs{param.binmode + 1};
    end
    if length(fms) ~= length(rate)
        annotation(gcf,'textbox','verticalalignment', 'bottom',...
            'String',{'Error in MTF'}, 'FontSize', 7, 'FitBoxToText','on');
        BMF = 100;
        WMF = 100;
        MTF_shape = 'Error';
        title(MTF_shape);
    else
        hold on
        line(ax, [1 fms(end)], [1 1]*rate(1),'Color',[0.7 0.7 0.7]); % noise alone line
        errorbar(ax, fms,rate, rate_std/sqrt(param.nrep),'.');
        line(ax, fms,rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w');
        plot(ax, fms,rate_sm,'-b', 'LineWidth', 1) % smoothed MTF

        % Null hypothesis: The spike rate with a tone at a specific
		j = NaN(num_mod_freqs, 1);
		p = NaN(num_mod_freqs, 1);
        for f_ind = 1:num_mod_freqs
            [j(f_ind),p(f_ind),~,~] = ttest(raw_rates(f_ind,1, :),raw_rates(1,1,:), 'Tail', 'both');
        end

        indices = find(p<0.05);
        scatter(fms(indices),rate(indices), 15, 'MarkerEdgeColor', [0.8500 0.3250 0.0980],'MarkerFaceColor', [0.8500 0.3250 0.0980]); % original results
        rat = squeeze(raw_rates(:, 1, :));
        for irep = 1:param.nrep
            scatter(fms(1:end,1), rat(:,irep), 3, 'k', 'filled');
        end
        hold off

        % Label the plots
        xtick = [1 2 5 10 20 50 100 200 500];
        xlim([fms(1) fms(end)])
        xlabel('Modulation Freq (Hz)')
        ylabel('Rate (sp/s)')
        set(ax,'XTick',xtick,'XScale', 'log', 'FontSize', 7)
        [BMF,WMF,MTF_shape, at_100] = MTFclassification(spike_rates, fms, fmi);

        if isnan(BMF)
            title(sprintf('%s, WMF = %dHz %s, No = %d',...
                MTF_shape, round(WMF), binmode, No));
        elseif isnan(WMF)
            title(sprintf('%s, BMF = %dHz %s, No = %d',...
                MTF_shape, round(BMF), binmode, No));
        elseif isnan(BMF) && isnan(WMF)
            title(sprintf('%s, %s, No = %d, No = , OnWin = %dms',...
                MTF_shape,  binmode, No, param.onsetWin));
        else
            title(sprintf('%s, WMF/BMF = %d/%dHz %s, No = %d',...
                MTF_shape, round(WMF),round(BMF), binmode, No));
        end

		% Create struct that contains processed data 
		data.rate = rate;
        data.fms = fms;
        data.rate_sm = rate_sm;
        data.rate_std = rate_std;
		data.t_spike = cluster.t_spike;
		data.t_spike_rel = cluster.t_spike_rel;
		data.stim_times = stim.times;
		data.abs_stim_num = cluster.abs_stim_num;
		data.at_100 = at_100;
    end
    grid on
    box on
    axis_set = axis;
    axis_set(3) = 0;
    axis(axis_set);
else
    MTF_shape = 'NaN';
    BMF = [];
    WMF = [];
    fig = [];
	data.at_100 = NaN;
end
data.MTF_shape = MTF_shape;
data.BMF = BMF;
data.WMF = WMF;

end