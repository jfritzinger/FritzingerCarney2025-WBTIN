function evaluate_fits(CF, data_rates, params, AN, model_params, ...
	fit_params_all)
%% evaluate_fits - Compare model predictions to experimental neural data
% 
% This function visualizes the performance of fitted IC models by comparing
% simulated responses to recorded spike rates across three stimulus types:
%   1. Modulation Transfer Function (MTF)
%   2. Narrowband Tuned Inhibitory Noise (TIN)
%   3. Wideband TIN (WB-TIN)
%
% Inputs:
%   CF            - Characteristic frequency of neuron (Hz)
%   data_rates    - Spike rates [MTF; TIN; WB-TIN]
%   params        - Cell array of stimulus parameters for each condition
%   AN            - Auditory nerve model outputs
%   model_params  - Lateral inhibition model parameters
%   fit_params_all- Matrix of fitted IC parameters (rows=CF ranges, cols=parameters)
%
% Outputs:
%   Generates tiled figures comparing data vs model predictions for all parameter sets
%
% Author: J. Fritzinger




%% Evaluate model fits  

% Plot data
figure('Position',[94,1,1404,946]);
tiledlayout(6, 6, 'padding', 'compact')
data_MTF = data_rates(1:26);
data_TIN = data_rates(27:29);
data_WB = data_rates(30:end);

num_paramCF = size(AN, 1);
for iparamCF = 1:num_paramCF

	% Run model with fit parameters
	AN_sub = AN(iparamCF,:);
	fit_params = fit_params_all(iparamCF,:);
	CS_params = [fit_params(1:2) 0];
	BMFs = fit_params(3:5);
	nstim = size(params, 2);
	model_outputs = cell(nstim, 1);
	for istim = 1:nstim
		param = params{istim};
		an_sout = squeeze(AN_sub{istim}.an_sout);
		an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
		an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
		model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
			an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
			'BMFs', BMFs);
	end

	% Analyze model output
	for istim = 1:nstim
		param = params{istim};
		model_output = model_outputs{istim};
		switch param.type
			case 'typMTFN'
				[~, model_MTF, ~, ~] = plotMTF(param, model_output.avIC, 0);
			case 'TIN'
				[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
			case 'SPEC_slide' % WB-TIN only for now
				[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
		end
	end

	% Plot MTF
	nexttile
	param = params{1};
	hold on
	yline(data_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,data_MTF);
	yline(model_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,model_MTF);
	%plot(data.fms,rate_sm,'-b', 'LineWidth', 1)
	hold off
	xtick = [1 2 5 20 50 100 200 500];
	xlim(xtick([1 end]))
	xlabel('Modulation Freq (Hz)')
	ylabel('Avg. Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	legend('Unmodulated', 'Location','best')
	grid on
	axis_set = axis;
	axis_set(3) = 0;
	axis(axis_set);
	title(sprintf('CF range: %d, [%0.02f %0.02f 0], BMF=[%0.0f %0.0f %0.0f]', ...
		iparamCF, fit_params(1), fit_params(2), fit_params(3),...
		fit_params(4), fit_params(5)))

	% Plot NB-TIN
	nexttile
	hold on
	num_SNRs = 3;
	SNRs = [-inf params{2}.SNR];
	x = repmat(1:num_SNRs, 1, 1);
	plot(x, data_TIN, 'LineWidth',1.5);
	plot(x,model_TIN, 'LineWidth',1.5);
	xlabel('SNR');
	xticks(1:num_SNRs)
	xticklabels([-inf 30 40])
	xlim([0 num_SNRs+1])
	ylabel('Spike rate (sp/s)')
	y = ylim;
	ylim([0,y(2)])
	legend('Data', 'Model', 'Location','best')

	% Plot WB-TIN
	nexttile
	param = params{3};
	freqs = param.fpeaks;
	hold on
	yline(mean(data_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,data_WB, 'LineWidth',2);
	yline(mean(model_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,model_WB, 'LineWidth',2);
	xline(CF/1000, '--', 'Color', [0.5 0.5 0.5]);
	hold off
	xlabel('Tone Frequency (kHz)')
	ylabel('Avg Rate (sp/s)')
	set(gca, 'XScale', 'log');
	grid on
	box on
	xlim([param.fpeaks(2) param.fpeaks(end)]/1000);
	axis_set = axis;
	axis_set(3) = 0;
	axis(axis_set);

end
