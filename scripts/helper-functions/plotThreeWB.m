function plotThreeWB(model, params, CF)
%PLOTMODELRESPONSES Visualizes auditory model responses across parameters
%   plotModelResponses(MODEL, PARAMS, CF, IWB) generates parameterized response
%   curves for auditory model analysis
%
%   Inputs:
%       model  - Nested cell array containing model data structures with:
%               .avIC : Average inferior colliculus response rates
%       params - Cell array of parameter structures containing:
%               .mlist    : Metadata with SNR/fpeak fields
%               .freq_lo  : Lower frequency bound (Hz)
%               .freq_hi  : Upper frequency bound (Hz)
%       CF     - Characteristic frequency (Hz) for reference line
%       iWB    - Index selector for parameter set (1-based)


linewidth = 1;
colors1 = {'#969696','#636363', '#252525'};
colors2 = {'#fd8d3c','#e6550d','#a63603'};


for iparam = 1:3
	lateral_model = model{iparam, 2};

	% Analyze model
	[SNRs,~,si] = unique([params{2}.mlist.SNR].');
	num_SNRs = length(SNRs);
	[fpeaks,~,fi] = unique([params{2}.mlist.fpeak].');
	num_fpeaks = length(fpeaks);
	rate_size = [num_fpeaks,num_SNRs];
	[avIC,~,~,~] = accumstats({fi,si},lateral_model.avIC, rate_size);

	% Plot model
	yline(avIC(1,1), 'Color', colors1{iparam}, 'LineWidth',linewidth)
	hold on
	plot(fpeaks/1000, avIC(:,2), 'LineWidth',linewidth, 'Color',colors2{iparam})

end
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
xlim([params{2}.freq_lo params{2}.freq_hi]./1000)
ylim([0,45])
set(gca, 'XScale', 'log');
grid on

end