function [DOGparams, R2, noise_alone_norm_all, max_rate_all] = ...
	fitDifferenceofGaussians(data, imsetype)


% DoG Fit
num_SNRs = length(data.SNRs);
R2 = zeros(num_SNRs);
noise_alone_norm_all = zeros(num_SNRs);
max_rate_all = zeros(num_SNRs);
DOGparams = zeros(num_SNRs, 6);
best_x = zeros(num_SNRs, 6);
for isnr = num_SNRs %1:num_SNRs

	rate = data.rate(2:end, isnr);
	noise_alone = mean(data.rate(1, isnr)); %:));
	r_fit = rate - noise_alone;

	% Normalize rates to a max of 1, rate array doesn't include noise_alone rate
	noise_alone = mean(data.rate(1, isnr)); %:));
	max_rate = max(data.rate(:,isnr));
	noise_alone_norm = noise_alone./max_rate;

	f_fit = data.fpeaks(2:end,1);
	f_fit = log10(f_fit);

	f_interp = linspace(f_fit(1), f_fit(end), 91)';
	r_fit2 = interp1(f_fit, r_fit,f_interp, 'linear');

	% FIX 
	tone_lo = data.fpeaks(2);
 	tone_hi = data.fpeaks(end);

	best_fval = Inf;
	for istarts = 1:100

		lb_CF = f_interp(1);
		ub_CF = f_interp(end);
		exc_s = 200 * rand(1);
		inh_s = 200 * rand(1);
		exc_init = 0.02 + (0.2 - 0.02) * rand(1);
		inh_init = 0.1 + (0.5 - 0.1) * rand(1);
		exc_CF_init = lb_CF + (ub_CF - lb_CF) * rand(1);
		inh_CF_init = lb_CF + (ub_CF - lb_CF) * rand(1);

		min_dif = f_fit(2)-f_fit(1);

		x0 = [exc_s, inh_s, exc_init, inh_init, exc_CF_init, inh_CF_init]; % Initial conditions
		lb = [0.1, 0.1, min_dif/2, min_dif/2, log10(tone_lo), log10(tone_lo)];
		ub = [200, 200, 3, 3, log10(tone_hi), log10(tone_hi)];

		% Parameter search
		fun = @(x) minimize_dog_fit(r_fit2, f_interp, x, imsetype);
		options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
			'MaxFunEvals', 10^12, 'maxiterations', 1500, 'ConstraintTolerance', 1e-12, ...
			'StepTolerance', 1e-16, 'display', 'off');
		[DOGparams(isnr,:), fval] = fmincon(fun, x0, [], [], [], [], lb, ub, [], options);

		if fval < best_fval
			best_x(isnr,:) = DOGparams(isnr,:);
			best_fval = fval;
		end
	end
	DOGparams = best_x;
	fprintf('Sexc=%0.2f, Sinh=%0.2f\n', DOGparams(isnr, 1), DOGparams(isnr, 2));


	% Calculate Coefficient of Determination
	[dog_rate, ~, ~] = createDoG(DOGparams(isnr,:), f_fit);
	dog_rate = (dog_rate+noise_alone_norm).*max_rate;
	R_int = corrcoef(dog_rate,data.rate(2:end, isnr));
	R2(isnr) = R_int(1, 2).^2;

	% Other saved data
	noise_alone_norm_all(isnr, :) = noise_alone_norm;
	max_rate_all(isnr) = max_rate;

end
end
