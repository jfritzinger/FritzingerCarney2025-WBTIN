function fig9_model_fit_BS(save_fig)
% model_fit_BS plots Figure 9.
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

% BS example 
putative = 'R29_TT4_P2_N2';
CF = 3482;

% Load in data and model
load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');
filename = sprintf('%s_ExpandedStimuli.mat', putative);
load(fullfile(datapath,putative, filename), 'SFIE', "broad_inh", "params")

%% Create figure

figure('Position',[13,603,4.567*ppi,3*ppi]);
legsize = 6;
titlesize = 8;
fontsize = 7;
spont_color = [0.4 0.4 0.4];
linewidth = 1;
labelsize = 13;

%% Plot example neuron

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
line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',linewidth);
errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.k', 'LineWidth',linewidth, 'CapSize',4);
line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
hold off
set(gca, 'XScale', 'log');
ylim([0 42])
xlim([data_MTF.fms(1) data_MTF.fms(end)])
xticks([1 2 5 10 20 50 100 200 500])
grid on
xticklabels([])
set(gca, 'fontsize', fontsize)
title('MTF', 'fontsize', titlesize);

% Plot NB-TIN
h(2) = subplot(4, 4, 2); % 11
rate_onCF = NaN(length(data_NB), 3);
rate_std_onCF = NaN(length(data_NB), 3);
label = cell(length(data_NB), 1);
for ind = 1:length(data_NB)
	NBTIN = data_NB{ind};
	param = params_NB{ind};
	[~,closestIndex] = min(abs(NBTIN.fpeaks(:, 1)-CF));
	SNRs = categorical(param.SNRNo);
	rate_onCF(ind, :) = NBTIN.rate(closestIndex,:);
	rate_std_onCF(ind,:) = NBTIN.rate_std(closestIndex,:);
	No = param.No;
	label{ind} = [num2str(No) ' dB SPL'];
end
errorbar(SNRs', rate_onCF', (rate_std_onCF/sqrt(param.nrep))', 'Linewidth', linewidth);
ylim([0 70])
xticklabels([])
grid on
set(gca, 'fontsize', fontsize)
title('NB-TIN','FontSize',titlesize);
box off

% Plot WB-TIN
iWB = 2;
h(3) = subplot(4, 4, 3); % 11
spont_color = [0.4 0.4 0.4];
colors = {'#82BB95', '#3F985C', '#034E1C'};
hold on
yline(mean(data_WB{iWB}.rate(1)),'Color',spont_color, 'LineWidth', linewidth);
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
errorbar(data_WB{iWB}.fpeaks(:,2),data_WB{iWB}.rate(:,2),data_WB{iWB}.rate_std(:,2)/(sqrt(params_WB{iWB}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{2});
set(gca, 'XScale', 'log');
xticklabels([])
ylim([0 80])
grid on
set(gca, 'fontsize', fontsize)
title('WB-TIN', 'fontsize', titlesize)

h(4) = subplot(4, 4, 4); % 11
colors = {'#82BB95', '#3F985C', '#034E1C'};
hold on
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.5); % CF line
errorbar(data_WB{1}.fpeaks(:,2),data_WB{1}.rate(:,2),data_WB{1}.rate_std(:,2)/(sqrt(params_WB{1}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{1});
errorbar(data_WB{2}.fpeaks(:,2),data_WB{2}.rate(:,2),data_WB{2}.rate_std(:,2)/(sqrt(params_WB{2}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{2});
errorbar(data_WB{3}.fpeaks(:,2),data_WB{3}.rate(:,2),data_WB{3}.rate_std(:,2)/(sqrt(params_WB{3}.nrep)), ...
	'LineWidth', linewidth,'Color', colors{3});
set(gca, 'XScale', 'log');
xticklabels([])
ylim([0 80])
grid on
set(gca, 'fontsize', fontsize)
title('All levels', 'FontSize',titlesize)

%% Plot lateral model

h(5) = subplot(4, 4, 5);
avIC = plotModelMTF(params{1}, broad_inh{1}.avIC);
ylabel('Avg. Rate (sp/s)')
set(gca, 'fontsize', fontsize)
ylim([0 40])
r = corrcoef(avIC,data_MTF.rate);
fprintf('Lateral Model, MTF R = %.2f\n', r(1,2))
grid on

h(6) = subplot(4, 4, 6);
avIC = plotModelTIN(params{2}, broad_inh{2}.avIC);
set(gca, 'fontsize', fontsize)
ylim([0 70])
r = corrcoef(avIC(2,:),rate_onCF);
fprintf('Lateral Model, NB-TIN R = %.2f\n', r(1,2))
grid on

h(7) = subplot(4, 4, 7);
avIC = plotModelWBTIN(params{4}, broad_inh{4}.avIC, CF, colors);
set(gca, 'fontsize', fontsize)
ylim([0 50])
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
ylim([0 50])
ylabel('Avg. Rate (sp/s)')
hLeg = legend('3 dB SPL', '23 dB SPL', '43 dB SPL', ...
	'Location','southwest', 'fontsize', legsize-2, 'box', 'off');
hLeg.ItemTokenSize = [6,6];
grid on


%% Plot SFIE
h(9) = subplot(4, 4, 9);
avIC = plotModelMTF(params{1}, SFIE{1}.average_ic_sout_BS);
ylim([0 40])
set(gca, 'fontsize', fontsize)
r = corrcoef(avIC,data_MTF.rate);
fprintf('SFIE, MTF R = %.2f\n', r(1,2))
xlabel('Mod. Freq. (Hz)')
xticks([1 2 5 10 20 50 100 200 500])
xticklabels([1 2 5 10 20 50 100 200 500])
hLeg = legend('Unmodulated', 'Location','southwest', 'fontsize', legsize,...
	'box', 'off');
hLeg.ItemTokenSize = [6,6];
grid on

h(10) = subplot(4, 4, 10);
avIC = plotModelTIN(params{2}, SFIE{2}.average_ic_sout_BS);
ylim([0 60])
r = corrcoef(avIC(2,:),rate_onCF);
fprintf('SFIE, NB-TIN R = %.2f\n', r(1,2))
xticklabels({'-Inf', '30', '40'})
xlabel('SNR')
set(gca, 'fontsize', fontsize)
hLeg = legend('23 dB SPL', 'Location','southwest', 'fontsize', legsize, ...
	'box', 'off');
hLeg.ItemTokenSize = [6,6];
grid on

h(11) = subplot(4, 4, 11);
avIC = plotModelWBTIN(params{4}, SFIE{4}.average_ic_sout_BS, CF, colors);
ylim([0 60])
set(gca, 'fontsize', fontsize)
r = corrcoef(avIC(:,2),data_WB{1}.rate(:,2));
fprintf('SFIE, WB-TIN, 3 dB, R = %.2f\n', r(1,2))
xlabel('Tone Freq. (kHz)')
xticks([0.5 1 2 3 4 6 8])
xticklabels([0.5 1 2 3 4 6 8])
hLeg = legend('Noise alone', '23 dB SPL', 'Location','southwest', ...
	'fontsize', legsize, 'box', 'off');
hLeg.ItemTokenSize = [6,6];
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
	[avIC,stdIC,~,~] = accumstats({fi,si},SFIE{2+ii}.average_ic_sout_BS, rate_size);
	errorbar(params{2+ii}.fpeaks/1000, avIC(:,2), stdIC(:,2)./sqrt(params{2+ii}.mnrep),...
		'LineWidth', linewidth, 'Color', colors{ii})
end
xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
ylim([0 60])
xlim([params{2+ii}.freq_lo params{2+ii}.freq_hi]./1000)
xticks([0.5 1 2 3 4 6 8])
xticklabels([0.5 1 2 3 4 6 8])
set(gca, 'XScale', 'log');
grid on
xlabel('Tone Freq. (kHz)')
set(gca, 'fontsize', fontsize)
grid on

%% Arrange Figure
height = 0.22;
width = 0.16;
left = [linspace(0.125, 0.52, 3) 0.815];
bottom = fliplr(linspace(0.11, 0.7, 3));

row = reshape(repmat(1:4, 4, 1), 16, 1);
col = repmat(1:4, 1, 4);
for ii = 1:12
	irow = row(ii);
	icol = col(ii);
	set(h(ii), 'position', [left(icol) bottom(irow) width height])
end

labels = {'Data', 'Broad Inh.', 'SFIE'};
labelbottom = [0.73 0.38 0.13];
for ii = 1:3
	annotation('textbox',[0.08 labelbottom(ii) 0.1 0.1],...
		'String',labels{ii},'FontWeight','bold','FontSize',titlesize,...
		'EdgeColor','none', 'Rotation', 90, 'FontName', 'Arial');
	annotation('textbox',[0.772 labelbottom(ii) 0.1 0.1],...
		'String',labels{ii},'FontWeight','bold','FontSize',titlesize,...
		'EdgeColor','none', 'Rotation', 90, 'FontName', 'Arial');
end

labels = {'A', 'B', 'C'};
labelbottom = fliplr(linspace(0.35, 0.96, 3));
for ii = 1:3
	annotation('textbox',[0.01 labelbottom(ii) 0.071 0.044],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end
annotation('textbox',[0.7 0.96 0.071 0.044],...
		'String','D','FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');

%% Save Figure

if save_fig == 1
	saveFigure('Fig9_model_fit_BS')
end
end
