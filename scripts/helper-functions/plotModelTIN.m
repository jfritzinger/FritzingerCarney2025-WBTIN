function avIC = plotModelTIN(params, IC_response)
%PLOTMODELTIN Visualizes TIN model responses
%   AVIC = PLOTMODELTIN(PARAMS, IC_RESPONSE) generates error bar plots of
%   model responses across different noise conditions and SNRs.
%
%   Inputs:
%       params - Struct containing stimulus parameters:
%           .mlist : 
%               .No  - Noise level
%               .SNR - Signal-to-Noise Ratio
%       IC_response - Model response data array
%
%   Output:
%       avIC - Averaged response matrix [NoiseConditions × SNRs]


linewidth = 1;

% Extract and categorize noise conditions and SNRs
[Nos,~,Noi] = unique(double([params.mlist.No]).'); % noise IDs for reproducible noises
num_Nos = length(Nos);
[SNRs,~,SNRi] = unique(double([params.mlist.SNR]).'); % noise IDs for reproducible noises
num_SNRs = length(SNRs);

% Calculate average responses and standard deviations
rate_size = [num_Nos,num_SNRs];
[avIC,rate_std] = accumstats({Noi,SNRi}, IC_response ,rate_size);

% Generate x-axis positions
x = repmat(1:num_SNRs, num_Nos, 1);

% Plot error bars 
if size(x, 1)>1
	errorbar(x(2,:)', avIC(2,:)', rate_std(2,:)', 'LineWidth',linewidth);
else
	errorbar(x', avIC', rate_std', 'LineWidth',linewidth);
end
xticks(1:num_SNRs)
xticklabels([])
xlim([0 num_SNRs+1])
colormap(jet)
box off
grid on
end