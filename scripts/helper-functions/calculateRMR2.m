function [R2_temp, BF_spl, iBF] = calculateRMR2(params_WB, data_WB, data_RM, snr_choice)

% Find BF nearest to WBTIN overall level
lbw = params_WB.fpeak_mid * 2^(-params_WB.bandwidth/2);
hbw = params_WB.fpeak_mid * 2^(params_WB.bandwidth/2);
overall_spl = params_WB.No + 10*log10(hbw-lbw);
[~,iBF] = min(abs(data_RM.SPL-overall_spl));
BF_spl = data_RM.SPL(iBF);

num_snr = length(snr_choice);
R2_temp = NaN(1, num_snr);
R2 = NaN(1, num_snr);
for i = 1:num_snr
	isnr = snr_choice(i);

	% Interpolate RM rate
	rate_WBTIN = data_WB.rate(2:end,isnr);
	fpeaks = data_WB.fpeaks(2:end,1);
	freqs_interp = logspace(log10(data_WB.fpeaks(2,1)), log10(data_WB.fpeaks(end,1)), 50);
	rate_interp_WBTIN = interp1(data_WB.fpeaks(2:end,1),rate_WBTIN,freqs_interp,'pchip'); % interpolates rate


	rate_RM = data_RM.rates(:,iBF);
	freqs_mid = logspace(log10(data_RM.freqs(1)), log10(data_RM.freqs(end)), 10000);
	rate_mid_RM = interp1(data_RM.freqs,rate_RM,freqs_mid,'pchip'); % interpolates rate

	% Uses interpolated WB-TIN
	[~, starting] = min(abs(freqs_mid-freqs_interp(1)));
	[~, ending] = min(abs(freqs_mid-freqs_interp(end)));
	rate_interp_RM = interp1(freqs_mid(starting:ending),rate_mid_RM(starting:ending),freqs_interp,'pchip'); % interpolates rate
	% Calculate correlation coefficient & variance explained
	R_int = corrcoef(rate_interp_RM,rate_interp_WBTIN);
	R2(i) = R_int(1, 2).^2;

	% Uses raw WB-TIN
	[~, starting] = min(abs(freqs_mid-fpeaks(1)));
	[~, ending] = min(abs(freqs_mid-fpeaks(end)));
	rate_RM_temp = interp1(freqs_mid(starting:ending),rate_mid_RM(starting:ending),fpeaks, 'pchip');
	
	% Calculate correlation coefficient & variance explained
	R_int = corrcoef(rate_RM_temp,rate_WBTIN);
	R2_temp(i) = R_int(1, 2).^2;

end
end