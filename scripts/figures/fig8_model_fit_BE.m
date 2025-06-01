function fig8_model_fit_BE(save_fig)
% model_fit_BE plots Figure 8.
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required:
%	 .mat files required:
%	 spreadsheets required: 
%
% Author: J. Fritzinger
% Created: 2025-01-20; Last revision: 2025-03-07
%
% -------------------------------------------------------------------------

%% Load in data 

[~, datapath, ~, ppi] = get_paths();

% BE example 
putative = 'R29_TT4_P2_N16';
CF = 2639;

% Load in data 
load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');

% Load in model 
filename = sprintf('%s_ExpandedStimuli.mat', putative);
load(fullfile(datapath, putative, filename),"rms_filt", 'SFIE', "broad_inh", "params")

%% Run & save models
% 
% addpath(fullfile(base, 'scripts', 'UR_EAR_2022a'), '-end');
% addpath(fullfile(base, 'scripts', 'UR_EAR_2022a', 'source'), '-end');
% addpath(fullfile(base, 'scripts', 'model-lat-inh'), '-end');
% addpath(fullfile(base, 'scripts', 'model-SFIE'), '-end');
% addpath(fullfile(base, 'scripts', 'model-energy'), '-end');
% 
% load(fullfile(datapath, putative, [putative '_BestModel.mat']), ...
% 	'fit_params', 'AN_best', 'params', 'model_params')
% 
% lm_params = [fit_params(1:2) 0];
% BMFs = fit_params(3:5);
% num_stim = size(params,2);
% for ist = 1:num_stim
% 	broad_inh{ist} = modelLateralSFIE_BMF(params{ist}, ...
% 		model_params, AN_best{ist}.an_sout, AN_best{ist}.an_sout_lo, ...
% 		AN_best{ist}.an_sout_hi,'CS_params', lm_params, 'BMF', BMFs);
% end
% 
% % Run SFIE Model
% model_params.CF_range = CF;
% model_params.CFs = CF;
% model_params.type = 'SFIE';
% for ist = 1:num_stim
% 	SFIE{ist} = wrapperIC(AN_best{ist}.an_sout, params{ist}, model_params);
% end
% 
% % Run Energy Model
% for ist = 1:num_stim
% 	param = params{ist};
% 	Fs = param.Fs;
% 	stimulus = [param.stim zeros(size(param.stim,1),0.1*Fs)];
% 	clear mdB;
% 	gamma_param.srate = Fs;
% 	tvals = (1:length(stimulus))/Fs;
% 	gamma_IF_reg = zeros(1,length(tvals));
% 	impaired = 0; % 0 = not impaired; 1 = 'impaired'
% 
% 	% Calculate RMS energy
% 	pin_gamma = zeros(size(stimulus, 1), Fs*param.dur+0.1*Fs);
% 	for istim = 1:size(stimulus, 1)
% 		gamma_param.fc = CF; % CF
% 		pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,...
% 			impaired, model_params.species);
% 
% 		% Analysis
% 		y2 = fft(stimulus(istim,:));
% 		m = abs(y2);
% 		mdB = 20*log10(m);
% 		f = (0:length(y2)-1)*Fs/length(y2);
% 		mdB(mdB<0) = 0;
% 		f(f>Fs/2) = [];
% 		mdB = mdB(1:length(f));
% 
% 	end
% 	pin_gamma = pin_gamma(:,1:param.dur*Fs);
% 	rms_filt{ist} = sqrt(mean(pin_gamma.^2,2));
% end
%
% % Save
% % save(fullfile(datapath, filename),...
% % 	'SFIE', "broad_inh", "rms_filt", "params")

%% Create figure

figure('Position',[50,50,4.567*ppi,3.8*ppi]);
legsize = 6;
titlesize = 8;
fontsize = 7;
spont_color = [0.4 0.4 0.4];
linewidth = 1;
labelsize = 13;

% Get data for each stimulus
params_WB = data(6:8,2); % Gets binaural WB-TIN stimuli
params_NB = data(10, 2);
params_MTF = data{3,2}; % Gets binaural MTFN
data_MTF = analyzeMTF(params_MTF);
data_WB = analyzeWBTIN(params_WB, []);
data_NB = analyzeNBTIN(params_NB, CF);

% Plot MTF
h(1) = subplot(4, 4, 1); %10
hold on
line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color,...
	'LineWidth',linewidth);
errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),...
	'.k', 'LineWidth',linewidth, 'CapSize',4);
line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',...
	5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
set(gca, 'XScale', 'log');
grid on
ylim([0 42])
xlim([data_MTF.fms(1) data_MTF.fms(end)])
xticks([1 2 5 10 20 50 100 200 500])
hLeg = legend('Unmodulated', 'Location','southwest', 'box','off');
hLeg.ItemTokenSize = [6,6];
xticklabels([])
set(gca, 'fontsize', fontsize)
title('MTF', 'fontsize', titlesize);

% Plot NB-TIN
h(2) = subplot(4, 4, 2); % 11
rate_onCF = NaN(length(data_NB), 3);
rate_std_onCF = NaN(length(data_NB), 3);
label = cell(length(data_NB), 1);
SNR = zeros(length(data_NB), 3);
for ind = 1:length(data_NB)
	NBTIN = data_NB{ind};
	param = params_NB{ind};
	[~,closestIndex] = min(abs(NBTIN.fpeaks(:, 1)-CF));
	SNR = categorical(param.SNRNo);
	rate_onCF(ind, :) = NBTIN.rate(closestIndex,:);
	rate_std_onCF(ind,:) = NBTIN.rate_std(closestIndex,:);
	No = param.No;
	label{ind} = [num2str(No) ' dB SPL'];
end
errorbar(SNR', rate_onCF', (rate_std_onCF/sqrt(param.nrep))', 'Linewidth', linewidth);
ylim([0 20])
xticklabels([])
set(gca, 'fontsize', fontsize)
grid on
title('NB-TIN','FontSize',titlesize);
box off

% Plot WB-TIN
iWB = 2;
h(3) = subplot(4, 4, 3); % 11
spont_color = [0.4 0.4 0.4];
colors = {'#82BB95', '#3F985C', '#034E1C'};
hold on
yline(mean(data_WB{iWB}.rate(1)),'Color',spont_color, 'LineWidth', linewidth);
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.5); % CF line
errorbar(data_WB{iWB}.fpeaks(:,2),data_WB{iWB}.rate(:,2),...
	data_WB{iWB}.rate_std(:,2)/(sqrt(params_WB{iWB}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{2});
set(gca, 'XScale', 'log');
xticklabels([])
ylim([0 42])
grid on
set(gca, 'fontsize', fontsize)
title('WB-TIN', 'fontsize', titlesize)

% Plot all levels WB-TIN
h(4) = subplot(4, 4, 4); % 11
colors = {'#82BB95', '#3F985C', '#034E1C'};
hold on
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.5); % CF line
errorbar(data_WB{1}.fpeaks(:,2),data_WB{1}.rate(:,2),...
	data_WB{1}.rate_std(:,2)/(sqrt(params_WB{1}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{1});
errorbar(data_WB{2}.fpeaks(:,2),data_WB{2}.rate(:,2),...
	data_WB{2}.rate_std(:,2)/(sqrt(params_WB{2}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{2});
errorbar(data_WB{3}.fpeaks(:,2),data_WB{3}.rate(:,2),...
	data_WB{3}.rate_std(:,2)/(sqrt(params_WB{3}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{3});
set(gca, 'XScale', 'log');
xticklabels([])
ylim([0 45])
grid on
set(gca, 'fontsize', fontsize)
title('All levels', 'fontsize', titlesize)

%% Plot lateral model

h(5) = subplot(4, 4, 5);
avIC = plotModelMTF(params{1}, broad_inh{1}.avIC);
ylabel('Avg. Rate (sp/s)')
set(gca, 'fontsize', fontsize)
ylim([0 15])
r = corrcoef(avIC,data_MTF.rate);
fprintf('Lateral Model, MTF R = %.2f\n', r(1,2))
grid on

h(6) = subplot(4, 4, 6);
avIC = plotModelTIN(params{2}, broad_inh{2}.avIC);
set(gca, 'fontsize', fontsize)
ylim([0 20])
r = corrcoef(avIC(2,:),rate_onCF);
fprintf('Lateral Model, NB-TIN R = %.2f\n', r(1,2))
grid on

h(7) = subplot(4, 4, 7);
avIC = plotModelWBTIN(params{4}, broad_inh{4}.avIC, CF, colors);
set(gca, 'fontsize', fontsize)
ylim([0 20])
r = corrcoef(avIC(:,2),data_WB{1}.rate(:,2));
fprintf('Lateral Model, WB-TIN, R = %.2f\n', r(1,2))
grid on

h(8) = subplot(4, 4, 8);
hold on
for ii = 1:3
	[SNRs,~,si] = unique([params{2+ii}.mlist.SNR].');
	num_SNRs = length(SNRs);
	[fpeaks,~,fi] = unique([params{2+ii}.mlist.fpeak].');
	num_fpeaks = length(fpeaks);
	rate_size = [num_fpeaks,num_SNRs];
	[avIC,stdIC,~,~] = accumstats({fi,si},broad_inh{2+ii}.avIC, rate_size);
	errorbar(params{2+ii}.fpeaks/1000, avIC(:,2), stdIC(:,2)./sqrt(params{2+ii}.mnrep),...
		'LineWidth', linewidth, 'Color', colors{ii})
end
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
xlim([params{2+ii}.freq_lo params{2+ii}.freq_hi]./1000)
xticklabels([])
set(gca, 'XScale', 'log');
grid on
set(gca, 'fontsize', fontsize)
ylim([0 25])
ylabel('Avg. Rate (sp/s)')


%% Plot SFIE
h(9) = subplot(4, 4, 9);
avIC = plotModelMTF(params{1}, SFIE{1}.average_ic_sout_BE);
ylim([0 40])
set(gca, 'fontsize', fontsize)
r = corrcoef(avIC,data_MTF.rate);
fprintf('SFIE, MTF R = %.2f\n', r(1,2))
grid on

h(10) = subplot(4, 4, 10);
avIC = plotModelTIN(params{2}, SFIE{2}.average_ic_sout_BE);
ylim([0 25])
set(gca, 'fontsize', fontsize)
r = corrcoef(avIC(2,:),rate_onCF);
fprintf('SFIE, NB-TIN R = %.2f\n', r(1,2))
grid on

h(11) = subplot(4, 4, 11);
avIC = plotModelWBTIN(params{4}, SFIE{4}.average_ic_sout_BE, CF, colors);
ylim([0 30])
set(gca, 'fontsize', fontsize)
r = corrcoef(avIC(:,2),data_WB{1}.rate(:,2));
fprintf('SFIE, WB-TIN, R = %.2f\n', r(1,2))
grid on

h(12) = subplot(4, 4, 12);
set(gca, 'fontsize', fontsize)
hold on
for ii = 1:3
	[SNRs,~,si] = unique([params{2+ii}.mlist.SNR].');
	num_SNRs = length(SNRs);
	[fpeaks,~,fi] = unique([params{2+ii}.mlist.fpeak].');
	num_fpeaks = length(fpeaks);
	rate_size = [num_fpeaks,num_SNRs];
	[avIC,stdIC,~,~] = accumstats({fi,si},SFIE{2+ii}.average_ic_sout_BE, rate_size);
	errorbar(params{2+ii}.fpeaks/1000, avIC(:,2), stdIC(:,2)./sqrt(params{2+ii}.mnrep),...
		'LineWidth', linewidth, 'Color', colors{ii})
end
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
ylim([0 44])
xlim([params{2+ii}.freq_lo params{2+ii}.freq_hi]./1000)
xticklabels([])
set(gca, 'XScale', 'log');
grid on
set(gca, 'fontsize', fontsize)

%% Plot Energy
h(13) = subplot(4, 4, 13);
avIC = plotModelMTF(params{1}, rms_filt{1});
xlabel('Mod. Freq. (Hz)')
xticks([1 2 5 10 20 50 100 200 500])
xticklabels([1 2 5 10 20 50 100 200 500])
ylim([0 0.1])
set(gca, 'fontsize', fontsize)
ylabel('RMS')
r = corrcoef(avIC,data_MTF.rate);
fprintf('Energy, MTF R = %.2f\n', r(1,2))
grid on

h(14) = subplot(4, 4, 14);
avIC = plotModelTIN(params{2}, rms_filt{2});
ylim([0 0.05])
xlabel('SNR');
xticklabels({'-Inf', '30', '40'})
set(gca, 'fontsize', fontsize)
hLeg = legend(label, 'Location', 'best', 'FontSize',legsize, 'box','off');
hLeg.ItemTokenSize = [6,6];
r = corrcoef(avIC(2,:),rate_onCF);
fprintf('Energy, NB-TIN R = %.2f\n', r(1,2))
grid on

h(15) = subplot(4, 4, 15);
avIC = plotModelWBTIN(params{4},rms_filt{4} , CF, colors);
ylim([0 0.05])
xticks([1 2 3 4])
xticklabels([1 2 3 4])
set(gca, 'fontsize', fontsize)
hLegend = legend({'Noise Alone', '23 dB SPL', 'CF'}, 'location','northwest',...
	'FontSize',legsize, 'box','off');
hLegend.ItemTokenSize = [6,6];
r = corrcoef(avIC(:,2),data_WB{1}.rate(:,2));
fprintf('Energy, WB-TIN, R = %.2f\n', r(1,2))
xlabel('Tone Freq. (kHz)')
grid on

h(16) = subplot(4, 4, 16);
set(gca, 'fontsize', fontsize)
hold on
for ii = 1:3
	[SNRs,~,si] = unique([params{2+ii}.mlist.SNR].');
	num_SNRs = length(SNRs);
	[fpeaks,~,fi] = unique([params{2+ii}.mlist.fpeak].');
	num_fpeaks = length(fpeaks);
	rate_size = [num_fpeaks,num_SNRs];
	[avIC,stdIC,~,~] = accumstats({fi,si},rms_filt{ii+2}, rate_size);
	errorbar(params{2+ii}.fpeaks/1000, avIC(:,2), stdIC(:,2)./sqrt(params{2+ii}.mnrep),...
		'LineWidth', linewidth, 'Color', colors{ii})
end
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
ylim([0 0.3])
xlim([params{2+ii}.freq_lo params{2+ii}.freq_hi]./1000)
set(gca, 'XScale', 'log');
grid on
ylabel('RMS')
hLegend = legend({'3 dB', '23 dB', '43 dB'}, 'location',...
	'northwest', 'FontSize',legsize-2, 'box','off');
hLegend.ItemTokenSize = [6,6];
xlabel('Tone Freq. (kHz)')
set(gca, 'fontsize', fontsize)

%% Arrange Figure
height = 0.16;
width = 0.16;
left = [linspace(0.13, 0.53, 3) 0.82];
bottom = fliplr(linspace(0.1, 0.78, 4));

row = reshape(repmat(1:4, 4, 1), 16, 1);
col = repmat(1:4, 1, 4);
for ii = 1:16
	irow = row(ii);
	icol = col(ii);
	set(h(ii), 'position', [left(icol) bottom(irow) width height])
end

labels = {'Data', 'Broad Inh.', 'SFIE', 'Energy'};
labelbottom = [0.76, 0.505 0.33 0.08];
for ii = 1:4
	annotation('textbox',[0.075 labelbottom(ii) 0.1 0.1],...
		'String',labels{ii},'FontWeight','bold','FontSize',titlesize,...
		'EdgeColor','none', 'Rotation', 90, 'FontName', 'Arial');
	annotation('textbox',[0.775 labelbottom(ii) 0.1 0.1],...
		'String',labels{ii},'FontWeight','bold','FontSize',titlesize,...
		'EdgeColor','none', 'Rotation', 90, 'FontName', 'Arial');
end

labels = {'A', 'B', 'C', 'D'};
labelbottom = fliplr(linspace(0.275, 0.97, 4));
for ii = 1:4
	annotation('textbox',[0.0000 labelbottom(ii) 0.071 0.044],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end
annotation('textbox',[0.72 0.96 0.071 0.044],...
		'String','E','FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');

%% Save Figure

if save_fig == 1
	saveFigure('Fig8_model_fit_BE')
end
end