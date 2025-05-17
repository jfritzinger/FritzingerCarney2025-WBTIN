function avIC = plotModelWBTIN(params, IC_response, CF, colors)
% Visualizes WBTIN model responses
%
%   Inputs:
%       params - Struct containing stimulus parameters:
%           .mlist : 
%               .SNR - Signal-to-Noise Ratio
%			.freq_lo
%			.freq_hi
%			.fpeaks
%       IC_response - Model response data array
%		CF: Neuron's CF
%		Colors: Cell array of colors for plot
%
%   Output:
%       avIC - Averaged response matrix [NoiseConditions × SNRs]


linewidth = 1;

% Analysis
[SNRs,~,si] = unique([params.mlist.SNR].');
num_SNRs = length(SNRs);
[fpeaks,~,fi] = unique([params.mlist.fpeak].');
num_fpeaks = length(fpeaks);
rate_size = [num_fpeaks,num_SNRs];
[avIC,~,~,~] = accumstats({fi,si},IC_response, rate_size);

% Plot
yline(avIC(1,1), 'k', 'LineWidth',1.5)
hold on
plot(params.fpeaks/1000, smooth(avIC(:,2)),...
	'LineWidth', linewidth, 'Color', colors{2})
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
xlim([params.freq_lo params.freq_hi]./1000)
xticklabels([])
set(gca, 'XScale', 'log');
end
