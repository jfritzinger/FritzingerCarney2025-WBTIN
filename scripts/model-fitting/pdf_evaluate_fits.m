function pdf_evaluate_fits(putative, CF, data_rates, params, AN, ...
	model_params, fit_params_all, modelpath)
% pdf_evaluate_fits - Generate a PDF report evaluating model fits to neural data
%
% Inputs:
%   putative      - (char) Neuron identifier string
%   CF            - (double) Characteristic frequency of the neuron (Hz)
%   data_rates    - (vector) Experimental spike rates for all conditions
%   params        - (cell array) Stimulus parameter structures
%   AN            - (cell array) Auditory nerve model outputs
%   model_params  - (struct) Parameters for the lateral inhibition model
%   fit_params_all- (matrix) Fitted IC model parameters (one row per CF range)
%   modelpath     - (char) Path to save the generated PDF report
%
% Outputs:
%   None (saves a PDF report to disk and displays it upon completion)
%
% Requirements:
%   - MATLAB Report Generator Toolbox (mlreportgen.dom, mlreportgen.report)
%   - Helper functions: modelLateralSFIE_BMF, plotMTF, plotTIN, plotWBTIN
%
% Author: J. Fritzinger


%% Set up PDF

% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
filename = 'Pearson';
report_name = fullfile(modelpath, sprintf('%s_%s.pdf', putative, filename));
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

% Evaluate model fits  
data_MTF = data_rates(1:26);
data_TIN = data_rates(27:29);
data_WB = data_rates(30:end);

num_paramCF = size(AN, 1);
for iparamCF = 1:num_paramCF

	% Run model with fit parameters
	AN_sub = AN(iparamCF,:);
	fit_params = fit_params_all(iparamCF,:);
	CS_params = [fit_params(1:2) 0.001];
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
	% iWB = 1;
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

	% Calculate MSE 
	mse = minimize_IC_model_fit(data_rates, AN_sub, params, model_params,...
		fit_params);
	% mse_all(iparamCF) = mse;

	% Paragraph 
	msg = sprintf(['CF range: %0.2f, [%0.02f %0.02f 0.001], ...' ...
		'BMF=[%0.0f %0.0f %0.0f], R=%0.2f'], ...
		AN_sub{1,1}.CF_span, fit_params(1), fit_params(2), fit_params(3),...
		fit_params(4), fit_params(5), mse);
	p = Paragraph(msg);
	p.FontSize = "14pt";
	p.WhiteSpace = "preserve";
	append(rpt,p);

	fig = figure();
	tiledlayout(1, 3, 'Padding','compact')

	% Plot MTF
	nexttile
	param = params{1};
	hold on
	yline(data_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,data_MTF);
	yline(model_MTF(1),'Color',[0.4 0.4 0.4], 'linewidth', 1.5);
	plot(param.all_fms,model_MTF);
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
	r = corrcoef(model_MTF, data_MTF);
	message = sprintf('R=%0.2f', r(1,2));
	text(0.05, 0.95, message, 'Units', 'normalized', ...
			'VerticalAlignment', 'top')

	% Plot NB-TIN
	nexttile
	hold on
	num_SNRs = 3;
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
	r = corrcoef(model_TIN, data_TIN);
	message = sprintf('R=%0.2f', r(1,2));
	text(0.05, 0.95, message, 'Units', 'normalized', ...
			'VerticalAlignment', 'top')

	nexttile
	param = params{3};
	freqs = param.fpeaks;
	hold on
	yline(mean(data_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,data_WB,  'LineWidth',2, 'Color',"#0072BD")
	yline(mean(model_WB(1)),'Color','k', 'LineWidth',2);
	plot(freqs/1000,model_WB,  'LineWidth',2, 'Color',"#D95319");
	xline(CF/1000, ':', 'Color', [0.5 0.5 0.5]);
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
	r = corrcoef(model_WB, data_WB);
	message = sprintf('R=%0.2f', r(1,2));
	text(0.05, 0.95, message, 'Units', 'normalized', ...
		'VerticalAlignment', 'top')

	tempname = ['CF' num2str(iparamCF)];
	[plt, images] = addToPDF(images, fig, tempname, [8, 2.1], modelpath);
	append(rpt, plt);

	p = Paragraph('  ');
	p.FontSize = "14pt";
	p.WhiteSpace = "preserve";
	append(rpt,p);

end

% Close report
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)


%% Functions

function [img, images] = addToPDF(images, fig, title, size, modelpath)
import mlreportgen.dom.*

% Set figure size, recommended
fig.PaperSize = size;
fig.PaperPosition = [0 0 size];
fig.Units = 'inches';
fig.Position(3:4) = size;

% Add the plot to the document
name = sprintf('%s.svg', title);
tempfilepath = fullfile(modelpath, name);
print(fig, tempfilepath, '-dsvg');
img = Image(tempfilepath);
delete(fig) %delete plot figure window
images = [images {img}];

end

end