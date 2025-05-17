function supp6_dog_analysis(save_fig)
% supp_dog_analysis plots Figure S6.
%	This function plots DoG analysis scatter plots for levels 3 and 43 dB
%	SPL. 
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN
%	 .mat files required: -
%	 spreadsheets required: DoG_Analysis.xlsx
%
% Author: J. Fritzinger
% Created: 2025-03-03; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------


%% Load in data

[~, datapath, ~, ppi] = get_paths();
analysisTable = readtable(fullfile(datapath, 'DoG_Analysis.xlsx'), ...
	'PreserveVariableNames',true);

%% Set up figure

figure('Name', 'DoG Analysis', 'Position',[100,100,6.9*ppi,6*ppi]);
legsize = 7;
titlesize = 9;
fontsize = 8;
labelsize = 13;
scattersize = 15;

%% Seperate data into groups

% Separate out 40 dB SNR
suprathreshold_ind = analysisTable.SNR==40;

% Seperate into 3, 23, and 43 No
No3 = analysisTable.No==0 | analysisTable.No==3;
No23 = analysisTable.No==20 | analysisTable.No==23;
No43 = analysisTable.No==40 | analysisTable.No==43;

% Good fit
good_fit = analysisTable.R2>=0.5;
good_unit = analysisTable.V_p>=0.4;

% Binmode
binaural = analysisTable.Binmode==2;

% BW = 4 octaves
noise_bw = analysisTable.NoiseBW==4 | analysisTable.NoiseBW==3;

%% Plot Fits

threshold_ind = analysisTable.SNR==40;
the = '40';

% 3 dB SPL
No = analysisTable.No==0 | analysisTable.No==3;
data_good_ind = threshold_ind & No & noise_bw & good_unit;
data_bad_ind = threshold_ind & No & noise_bw & ~good_unit;
R2_values = analysisTable.R2(data_good_ind);
R2_values_bad = analysisTable.R2(data_bad_ind);
h(1) = subplot(4, 4, 1);
edges = linspace(0, 1, 21);
hold on
histogram(R2_values_bad, edges,'FaceColor', "#D95319")
histogram(R2_values, edges, 'FaceColor',"#0072BD");
xlabel('R^2')
set(gca, 'FontSize', fontsize)
title([the ' dB SNR'], 'fontsize', titlesize)
xticks(0:0.2:1)
xlim([0 1])
hold on
ylim([0 80])
ylabel('# Fits')
hLegend = legend('V_p < 0.4','V_p > 0.4', 'location', 'northwest', 'fontsize', legsize);
hLegend.ItemTokenSize = [6,6];
title('N_0 = 3 dB SPL', 'FontSize',titlesize)

% 43 dB SPL
No = analysisTable.No==40 | analysisTable.No==43;
data_good_ind = threshold_ind & No & noise_bw & good_unit;
data_bad_ind = threshold_ind & No & noise_bw & ~good_unit;
R2_values = analysisTable.R2(data_good_ind);
R2_values_bad = analysisTable.R2(data_bad_ind);
h(2) = subplot(4, 4, 2);
edges = linspace(0, 1, 21);
hold on
histogram(R2_values_bad, edges,'FaceColor', "#D95319")
histogram(R2_values, edges, 'FaceColor',"#0072BD");
xlabel('R^2')
set(gca, 'FontSize', fontsize)
title([the ' dB SNR'], 'fontsize', titlesize)
xticks(0:0.2:1)
xlim([0 1])
hold on
ylim([0 80])
yticklabels([])
hLegend = legend('V_p < 0.4','V_p > 0.4', 'location', 'northwest', 'fontsize', legsize);
hLegend.ItemTokenSize = [6,6];
title('N_0 = 43 dB SPL', 'FontSize',titlesize)


%% Plot Sigmas

No_names = {'3', '23', '43'};
No = [No3, No23, No43];
iii = 1;
sigma_exc = zeros(3, 2); sigma_inh = zeros(3,2);
index = [5 0 6; 9 0 10; 13 0 14];
for itype = 1:3
	binmode = binaural;
	for iNo = [1 3]
		No_name = No_names{iNo};

		toAnalyze = suprathreshold_ind & No(:, iNo) & good_fit & binmode & noise_bw & good_unit;
		ind = find(toAnalyze);
		CFs = analysisTable.CF(ind);
		if itype == 1
			sig = 10.^analysisTable.CF_exc(ind).*(10.^analysisTable.sigma_exc(ind)-1);
			sigma_exc(iNo,1:2) = [mean(sig./CFs) median(sig./CFs)];
		elseif itype == 2
			sig = 10.^analysisTable.CF_inh(ind).*(10.^analysisTable.sigma_inh(ind)-1);
			sigma_inh(iNo,1:2) = [mean(sig./CFs) median(sig./CFs)];
		else
			sig_exc = 10.^analysisTable.CF_exc(ind).*(10.^analysisTable.sigma_exc(ind)-1);
			sig_inh = 10.^analysisTable.CF_inh(ind).*(10.^analysisTable.sigma_inh(ind)-1);
			sig = sig_inh./sig_exc;

			fprintf('No=%s: %d/%d neurons above 1\n', No_names{iNo}, sum(sig>=1), length(sig))
		end

		% Fit linear regression line
		mdl = fitlm(log10(CFs), log10(sig));
		x = 0.3:0.5:10000;
		p(1) = mdl.Coefficients.Estimate(2,1);
		p(2) = mdl.Coefficients.Estimate(1,1);
		p(3) = mdl.Coefficients.pValue(2);
		p(4) = mdl.Rsquared.Ordinary;
		if itype == 1
			mdlfit2(iNo,:) = 10.^(p(1)*log10(x)+p(2));
		end
		mdlfit = 10.^(p(1)*log10(x)+p(2));

		h(index(itype, iNo)) = subplot(4, 4, index(itype, iNo));
		hold on
		if itype == 1
			if p(3) < 0.05
				plot(x/1000, mdlfit, 'r');
			else
				plot(x/1000, mdlfit, '--r');
			end
			scatter(CFs/1000, sig,scattersize, 'filled', 'r', 'MarkerEdgeColor', 'k')
		elseif itype == 2
			plot(x/1000, mdlfit2(iNo,:), 'r');
			if p(3) < 0.05
				plot(x/1000, mdlfit, 'b');
			else
				plot(x/1000, mdlfit, '--b');
			end
			scatter(CFs/1000, sig,scattersize, 'filled', 'b', 'MarkerEdgeColor', 'k')
		else
			if p(3) < 0.05
				plot(x/1000, mdlfit, 'color', '#009E73');
			else
				plot(x/1000, mdlfit, '--', 'color', '#009E73');
			end
			scatter(CFs/1000, sig,scattersize, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', '#009E73')
		end

		if itype == 3
			ylim([0.1 10])
		else
			ylim([40 15000])
		end

		hold on
		xlim([0.3 11])
		set(gca, 'FontSize', fontsize)
		set(gca, 'yscale', 'log')
		set(gca, 'xscale', 'log')
		if itype == 2
			hLegend = legend('', '', sprintf('p = %.3f', p(3)), ...
				'Location', 'NorthWest', 'fontsize', legsize);
			hLegend.ItemTokenSize = [6,6];
		else
			hLegend = legend('', sprintf('p = %.3f', p(3)), ...
				'Location', 'NorthWest', 'fontsize', legsize);
			hLegend.ItemTokenSize = [6,6];
		end
		xticks([100 200 500 1000 2000 5000 10000]./1000)

		grid on

		if itype == 3
			xlabel('CF (kHz)')
		end

		if itype == 1
			title(sprintf('N_0 = %s dB SPL', No_name), 'FontSize',titlesize)
			xticklabels([])
			yticks([100 200 500 1000 2000 5000 10000 20000])
			yticklabels([100 200 500 1000 2000 5000 10000 20000]./1000)
		elseif itype == 2
			xticklabels([])
			yticks([100 200 500 1000 2000 5000 10000 20000])
			yticklabels([100 200 500 1000 2000 5000 10000 20000]./1000)
		elseif itype == 3
			yticks([0.1 0.2 0.5 1 2 5])
		end

		if iNo == 1
			if itype == 1
				ylabel('\sigma_e_x_c (kHz)')
			elseif itype == 2
				ylabel('\sigma_i_n_h (kHz)')
			else
				ylabel('\sigma_i_n_h/\sigma_e_x_c')
			end
		else
			yticklabels([])
		end

		iii = iii+1;
	end
end


%% CF Scatter plots

binmode = binaural;
index = [7 0 8; 11 0 12; 15 0 16];
for itype = 1:3
	for iNo = [1 3]
		toAnalyze = suprathreshold_ind & No(:,iNo) & good_fit & binmode & noise_bw;
		ind = find(toAnalyze);
		CFs = analysisTable.CF(ind);

		No_name = No_names{iNo};
		if itype == 1
			CF_exc = 10.^analysisTable.CF_exc(ind);
		elseif itype == 2
			CF_inh = 10.^analysisTable.CF_inh(ind);
		else
			CF_exc = 10.^analysisTable.CF_exc(ind);
			CF_inh = 10.^analysisTable.CF_inh(ind);
			CF_offset = log2(CF_exc) - log2(CF_inh); %(CF_exc - CF_inh)./analysisTable.CF(ind);
		end

		h(index(itype, iNo)) = subplot(4, 4, index(itype, iNo));
		hold on
		if itype == 1
			plot([min(CFs) max(CFs)]/1000, [min(CFs) max(CFs)]/1000, 'k')
			scatter(CFs/1000, CF_exc/1000, scattersize, 'filled', 'r', 'MarkerEdgeColor', 'k');
		elseif itype == 2
			plot([min(CFs) max(CFs)]/1000, [min(CFs) max(CFs)]/1000, 'k')
			scatter(CFs/1000, CF_inh/1000, scattersize, 'filled', 'b','MarkerEdgeColor', 'k');
		elseif itype == 3
			yline(0, 'k')
			scatter(CFs/1000, CF_offset,scattersize, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', '#009E73');
		end

		set(gca, 'FontSize', fontsize)
		xticks([10 20 50 100 200 500 1000 2000 5000 10000]./1000)
		set(gca, 'xscale', 'log')
		xlim([0.3 11])

		if itype == 1
			title(sprintf('N_0 = %s dB SPL', No_name), 'FontSize',titlesize)
			yticks([10 20 50 100 200 500 1000 2000 5000 10000]./1000)
			ylim([0.3 17])
			set(gca, 'yscale', 'log')
			xticklabels([])
		elseif itype == 2
			yticks([10 20 50 100 200 500 1000 2000 5000 10000]./1000)
			ylim([0.3 17])
			set(gca, 'yscale', 'log')
			xticklabels([])
		elseif itype == 3
			ylim([-1 1])
			xlabel('CF (kHz)')
		end
		grid on

		if iNo == 1
			if itype == 1
				ylabel('{\it f}_{exc} (kHz)')
			elseif itype == 2
				ylabel('{\it f}_{inh} (kHz)')
			elseif itype == 3
				ylabel('{\it f}_{exc} - {\it f}_{inh}')
			end
		else
			yticklabels([])
		end
	end
end

%% Strengths

No_names = {'3', '23', '43'};
No = [No3, No23, No43];

iii = 1;
binmode = binaural;
for iNo = [1 3]
	No_name = No_names{iNo};
	toAnalyze = suprathreshold_ind & No(:, iNo) & good_fit & binmode & noise_bw & good_unit;
	ind = find(toAnalyze);
	sig = analysisTable.g_ratio(ind);

	% Histogram of g_ratio
	if iNo == 1
		h(3) = subplot(4, 4, 3);
	elseif iNo == 3
		h(4) = subplot(4, 4, 4);
	end
	edges = linspace(0, 1.6, 25);
	histogram(sig, edges)
	hold on
	xline(1, 'k', 'linewidth', 1.5)
	xlim([0 1.7])
	ylim([0 20])
	if iNo == 1
		ylabel('# neurons')
	else
		yticklabels([])
	end
	xlabel('Strength Ratio (Inh/Exc)')
	set(gca, 'FontSize', fontsize)
	title(sprintf('N_0 = %s dB SPL', No_name), 'FontSize',titlesize)
	box off
	iii = iii+1;
end

%% Annotations / Rearrange figure

cols = [0.08 0.29 0.58 0.79];
heights = [0.06 0.27 0.48 0.81];
width = 0.2;
height = 0.2;

set(h(1), 'Position', [cols(1) heights(4) width height-0.06]) % [left, bottom, width, height]
set(h(2), 'Position', [cols(2) heights(4) width height-0.06]) % [left, bottom, width, height]

set(h(3), 'Position', [cols(3) heights(4) width height-0.06])
set(h(4), 'Position', [cols(4) heights(4) width height-0.06])

set(h(5), 'Position', [cols(1) heights(3) width height])
set(h(6), 'Position', [cols(2) heights(3) width height])
set(h(7), 'Position', [cols(3) heights(3) width height])
set(h(8), 'Position', [cols(4) heights(3) width height])

set(h(9), 'Position', [cols(1) heights(2) width height])
set(h(10), 'Position', [cols(2) heights(2) width height])
set(h(11), 'Position', [cols(3) heights(2) width height])
set(h(12), 'Position', [cols(4) heights(2) width height])

set(h(13), 'Position', [cols(1) heights(1) width height])
set(h(14), 'Position', [cols(2) heights(1) width height])
set(h(15), 'Position', [cols(3) heights(1) width height])
set(h(16), 'Position', [cols(4) heights(1) width height])

% Create textbox
annotation('textbox',[0.02 0.97 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.51 0.97 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.02 0.71 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.51 0.71 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');


%% Export Figure

if save_fig == 1
	saveFigure('Supp6_DoGAnalysis')
end
end