function [fig, rate,rate_std] = plotTIN(stim_params, model, plot_on)


[Nos,~,Noi] = unique(double([stim_params.mlist.No]).'); % noise IDs for reproducible noises
%Noi = repmat(Noi, stim_params.mnrep, 1);

%[SPLs,~,SPLi] = unique(double([stim_params.mlist.noise_spl]).');
num_Nos = length(Nos);
[SNRs,~,SNRi] = unique(double([stim_params.mlist.SNR]).'); % noise IDs for reproducible noises
%SNRi = repmat(SNRi, stim_params.mnrep, 1);

num_SNRs = length(SNRs);

rate_size = [num_Nos,num_SNRs];
[rate,rate_std] = accumstats({Noi,SNRi}, model ,rate_size);

% Plot
if plot_on == 1
	fig = figure;
	x = repmat(1:num_SNRs, num_Nos, 1);
	errorbar(x', rate', rate_std', 'LineWidth',1.5);
	for ii = 1:num_Nos
		%label{ii} = ['No = ' num2str(Nos(ii)) ' dB SPL, Overall level = '  num2str(round(SPLs(ii)))];
		label{ii} = ['No = ' num2str(Nos(ii))];
	end
	%hLegend = legend(label, 'Location', 'Southwest', 'FontSize',8);
	%hLegend.ItemTokenSize = [6,6];
	xlabel('SNR');
	xticks(1:num_SNRs)
	xticklabels(strtrim(num2str(SNRs)))
	xlim([0 num_SNRs+1])
	ylabel('Spike rate (sp/s) +/1 STD')
	y = ylim;
	ylim([0,y(2)])

	discrp = [];
	if stim_params.choice == 1
		discrp = 'Frozen noise, ';
	elseif stim_params.choice == 2
		discrp = 'Warm noise, ';
	elseif stim_params.choice == 3
		discrp = 'All warm noise, ';
	end

	if stim_params.binmode == 1
		discrp = [discrp 'Contra, '];
	elseif stim_params.binmode == 2 && stim_params.condition == 0
		discrp = [discrp 'NoSo, '];
	elseif stim_params.binmode == 2 && stim_params.condition == 1
		discrp = [discrp 'NoSpi, '];
	end

	if stim_params.bw_choice == 1
		discrp = [discrp 'NB, '];
	elseif stim_params.bw_choice == 2
		discrp = [discrp 'WB, '];
	elseif stim_params.bw_choice == 3
		discrp = [discrp '100Hz, '];
	elseif stim_params.bw_choice == 4
		discrp = [discrp '3 Oct, '];
	end

	title(['Tone frequency is ' num2str(stim_params.CF)...
		' Hz, ' discrp num2str(stim_params.mnrep) ' reps' ]);
	colormap(jet)
else
	fig = 0;
end

end