function SNR_er = snr_er_est(y)
% function from Pospisil and Bair 2021, python code translated into MATLAB
% Approximately unbiased estimator of snr.
%
%     Assumes y has equal variance across trials and observations
%     (may require variance stabilizing transform).
%
% Parameters
% ----------
% y : numpy.ndarray
%     N neurons X n trials X m observations array
%
% Returns
% -------
% snr_ER : an approximately unbiased estimate of snr
% --------

[n, m] = size(y);

% estimate of trial-to-trial variability
sig2_hat = mean(var(y));

% Approx unbiased estimator of dynamic range across expected values.
ym = mean(y);
y_ms = ym - mean(ym);
d2 = sum(y_ms.^2);
d2_er = d2 - (m-1)*sig2_hat/n;

%     if scale:
%         d2_er = d2_er/m


SNR_er = d2_er/sig2_hat;

end



