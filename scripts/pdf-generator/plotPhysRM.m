function [fig, data] = plotPhysRM(cluster, params, stims, CF)
% plotPhysRM plots the response map. 
%   Response map includes the average rate responses for all levels, and
%   the 'PortFors' map which shows excitation over the time of the stimulus
%   for all frequencies. 
%
% Inputs:
%    cluster - struct, output of post_process
%    params - params from post_process that have the type 'type=RM'
%	 stim - struct, output of post_process
%	 CF - CF of the neuron
%
% Outputs:
%    fig - Average rate plot
%    data - struct that contains processed average rate data. 
%		 data.rates - average rate of neuron
%        data.spont - spontaneous rate of neuron
%        data.freqs - tone frequencies from stimulus
%		 data.t_spike - spike times 
%	     data.stim_times - stimulus times 
%	     data.CF - CF of the neuron
%
% Other m-files required: accumstats, checkIfSPLShift
%
% Author: J. Fritzinger, adapted from cell_characterization
% Created: 2022-06-18; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------

% Process cluster 
if iscell(params)
	param = params{1}; % If multiple, plots only the first instance
else
	param = params;
end
ds = param.dsid;
dsi = 1;
if length(stims)==1
	this_ds = stims.dsid == ds(dsi);
	stim = stims;
else
	this_ds = stims(ds).dsid == ds(dsi);
	stim = stims(ds);
end


binmode_strs = {'Ipsilateral','Contralateral','Binaural'};
binmode = binmode_strs{param.binmode + 1};

[spls,~,si] = unique(double([param.list.spl]).');
num_spls = length(spls);
[freqs,~,fi] = unique(double([param.list.freq]).');
num_freqs = length(freqs);

% Check if this session has an spl shift
spl_shift = checkIfSPLShift(param.animal, param.session);
spls = spls + spl_shift;

rate_size = [num_freqs,num_spls];
num_spikes = cluster.num_spikes_peri(this_ds);
spike_rates = num_spikes*1000/param.dur; % spikes/sec

fig = figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
if length(si) == length(spike_rates)
    [rates,~,~,~] = accumstats({fi,si},spike_rates, rate_size);
    
    spont = mean(rates(:,1));       
    max_y = max(rates(:));
    if max_y == 0
        max_y = 5;
    end
    num_plot_rows = max(num_spls,5); % always plot at least 5 rows
    for spl_index = 1:num_spls
        subplot(num_plot_rows,2,2*num_plot_rows+1-spl_index*2)
        plot(log10(freqs([1 end])),[1 1]*spont,'-',...
            'Color',[0.75 0.75 0.75])
        set(gca,'FontSize',8)
        hold all

        if ~isnan(CF)
            xline(log10(CF), '--', 'Color', [0.5 0.5 0.5]);
        end
        
        plot(log10(freqs),rates(:,spl_index),'r-','LineWidth',1)   
            
        if spls(spl_index) == -Inf
            title(sprintf('SPL = %c%c dB (0 Pa), %s',8722,8734, binmode)) % minus infinity
        else
            title(sprintf('SPL = %.0f dB, %s',spls(spl_index), binmode))
        end
        set(gca,'XTick',log10([250 500 1e3 2e3 5e3 10e3]),...
            'XTickLabel',{'0.25','0.5','1','2','5','10'})
        xlims = log10(freqs([1 end]));
        xlim(xlims.*[0.98;1.02])
        ylim([0 max_y])
        text(xlims(2),0,' Freq. (kHz)','FontSize',8,'HorizontalAlignment','left',...
            'VerticalAlignment','top')
        ylabel({'rate','(sp/s)'})
        grid on

		% Create struct that contains processed data 
		data.SPL(spl_index) = spls(spl_index);
        [~, i_BF] = max(rates(:,spl_index));
        data.BF(spl_index) = freqs(i_BF);
        data.rates = rates;
        data.spont = spont;
        data.freqs = freqs;

    end
    hold off
    
    % Plot Portfors
    france = makecmap([0 0 1;1 1 1;1 0 0],256,'rgb');
    scale = -inf;
    
    tinc = 12500; % µs
    tinc_s = tinc/1e6; % seconds
    t_bins = 0:tinc:500000;
    num_bins = length(t_bins) - 1;
    portfors_bins = reshape(bsxfun(@plus,stim.times,t_bins).',[],1);
    
    [freqs,~,fi] = unique(double([param.list.freq]).');
    num_freqs = length(freqs);
    
    % Determine number of spikes per 12.5 ms bin.
    num_portfors_spikes = histc(cluster.t_spike,portfors_bins);
    num_portfors_spikes(41:41:end) = [];
    num_portfors_spikes = reshape(num_portfors_spikes,num_bins,[]);
    
    ti = 1:num_bins;
    [ti2,fi2] = ndgrid(ti,fi);
    siz = [num_bins,num_freqs];
    
    ax_handles = zeros(1,num_spls);
    for spl_index = 1:num_spls%:-1:1
        k = si == spl_index;
        ti3 = ti2(:,k);
        fi3 = fi2(:,k);
        N3 = num_portfors_spikes(:,k);
        Tc = accumarray({ti3(:),fi3(:)},tinc_s,siz);
        A = accumarray({ti3(:),fi3(:)},N3(:),siz)./Tc - spont;
        Asm = conv2(A,ones(3,3)/9,'same');
        
        ax_handles(spl_index) = subplot(num_plot_rows,6,...
            6*num_plot_rows+4-spl_index*6+[0 1]);
        pcolor(log10(freqs),t_bins(1:end-1)/1000,Asm)
        scale = max(scale,max(abs(Asm(:))));
        shading flat
        set(gca,'XTick',log10([250 500 1e3 2e3 5e3 10e3]),...
            'XTickLabel',{'0.25','0.5','1','2','5','10'},...
            'YDir','reverse','FontSize',8)
        ylabel({'time','(ms)'})
        colormap(france)
        colorbar
    end
    if scale == 0
        scale = 5;
    end
    set(ax_handles,'CLim',[-1 1]*scale)

	% Create struct that contains processed data 
	data.t_spike = cluster.t_spike;
	data.stim_times = stim.times;
	data.CF = CF;

else
    data = [];
end

end