function [fig, rate, rate_std] = plotWBTIN(stim_params, model_response, plot)

[SNRs,~,si] = unique([stim_params.mlist.SNR].');
num_SNRs = length(SNRs);

[fpeaks,~,fi] = unique([stim_params.mlist.fpeak].');
num_fpeaks = length(fpeaks);

rate_size = [num_fpeaks,num_SNRs];
[rate,rate_std,~,~] = accumstats({fi,si},model_response, rate_size);
freqs = repmat(fpeaks, 1, length(SNRs));

if plot == 1
	fig = figure;
	% BE Model response
	hold on
	line([freqs(2) freqs(end)]/1000, [1 1]*mean(rate(1)),'Color','k', 'LineWidth',2);
	errorbar(freqs/1000,rate,rate_std, 'LineWidth',2);
	xline(stim_params.CF/1000, '--', 'Color', [0.5 0.5 0.5]);
	%legend(char('Noise Alone', num2str(SNRs)))
	hold off
	xlabel('Tone Frequency (kHz)')
	ylabel('Avg Rate (sp/s)')
	set(gca, 'XScale', 'log');
	%ylim([0 70])
	grid on
	box on
	xlim([stim_params.fpeaks(2) stim_params.fpeaks(end)]/1000);
	axis_set = axis;
	axis_set(3) = 0;
	axis(axis_set);

	for ii = 1:num_SNRs
		leg(ii,:) = [num2str(round(SNRs(ii),1)) 'dB SNR'];
	end
	hLegend = legend(char('Noise Alone', leg, 'CF'), 'location', 'best', 'Fontsize', 8);
	hLegend.ItemTokenSize = [6,6];
else
	fig = [];
end
end