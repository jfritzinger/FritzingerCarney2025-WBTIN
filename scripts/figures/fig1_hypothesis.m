function fig1_hypothesis(save_fig)
% supp_dog_analysis plots Figure 1.
%	This function plots the model hypothesis for NB-TIN and WB-TIN
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN
%	 .mat files required: modelAN, wrapperIC, plotModelMTF, plotModelTIN, 
%							plotModelWBTIN
%	 spreadsheets required: -
%
% Author: J. Fritzinger
% Created: 2025-01-20; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------

%% Set up figure 

[~, ~, ~, ppi] = get_paths();
figure('Position',[50,50,4.567*ppi,2.2*ppi])
legsize = 6;
titlesize = 8;
fontsize = 7;
linewidth = 1;
labelsize = 13;

%% NB-TIN  

% Parameters
CF = 1200;
params.Fs = 100000;
params.mnrep = 1;
params.type = 'NB_noise';
params.version = 5;
params.fpeak_mid = CF;
params.ramp_dur = 0.01;
params.dur = 0.3;
params.reptim = 0.6;
params.CF = CF;
params.stp_otc = 0;
params.bandwidth = 0.3;
params.No = 40;
params.SNRNo = 40;
params.range = 2;
params.nrep = 1;
params.binmode = 2;
params.onsetWin = 25;
params.physio = 0;
params = generate_NBTIN(params);

% Analysis
params.stp_oct = 2;
params.freq_lo = params.fpeak_mid*2^(-params.range/2); % center-noise freqs
params.freq_hi = params.fpeak_mid*2^(params.range/2);

% Plot CF 
h(1) = subplot(3, 3, 1);
hold on
line([CF CF], [0 100], 'LineWidth', linewidth,'LineStyle', ':', 'Color',...
	[.7 .7 .7]);

% Plot on-CF 
lbw = params.CF * 2^(-params.bandwidth/2);
hbw = params.CF * 2^(params.bandwidth/2);
line([lbw hbw], [35 35], 'LineWidth', linewidth,'Color', '#014a2f')
line([CF CF], [35 75], 'LineWidth', linewidth,'Color', '#014a2f')
line([lbw lbw], [0 35.5], 'LineWidth', linewidth,'Color', '#014a2f')
line([hbw hbw], [0 35.5], 'LineWidth', linewidth,'Color', '#014a2f')

% Figure Properties
ylimits = [0 90];
xlimits = [275 4900];
set(gca, 'XScale', 'log')
ylim(ylimits)
yticks([])
xlim(xlimits)
xticks([600 1200 2400])
xticklabels([])
hleg = legend('CF', 'fontsize', legsize, 'box', 'off');
hleg.ItemTokenSize = [10, 8];
ylabel('Mag.')
set(gca,'fontsize',fontsize)
title('NB-TIN Stimulus', 'FontSize',titlesize)

%% WB-TIN

% Parameters
CF = 1200;
params.Fs = 100000;
params.mnrep = 1;
params.type = 'WB_noise';
params.version = 5;
params.fpeak_mid = CF;
params.ramp_dur = 0.01;
params.dur = 0.3;
params.reptim = 0.6;
params.CF = CF;
params.stp_oct = 0;
params.bandwidth = 4;
params.No = 40;
params.SNR = 40;
params.range = 3;
params.nrep = 1;
params.binmode = 2;
params.onsetWin = 25;
params.physio = 0;
params = generate_WBTIN(params);

% Analysis 
params.stp_oct = 6;
params.freq_lo = params.fpeak_mid*2^(-params.range/2); % center-noise freqs
params.freq_hi = params.fpeak_mid*2^(params.range/2); 
nfreqs = floor(params.stp_oct * log2(params.freq_hi/params.freq_lo))+1; 
fpeaks = [0 params.freq_lo * 2.^((0:nfreqs-1)/params.stp_oct)];
lbw = params.CF * 2^(-params.bandwidth/2);
hbw = params.CF * 2^(params.bandwidth/2);

% Plotting
h(2) = subplot(3, 3, 2);
hold on
line([CF CF], [0 100], 'LineWidth', linewidth,'LineStyle', ':', ...
	'Color', [.7 .7 .7]);
line([0.1  300], [35 35], 'LineWidth', linewidth,'LineStyle', ':',...
	'Color', [.7 .7 .7]);
for ind = 1:nfreqs
    line([fpeaks(ind+1) fpeaks(ind+1)], [35 75], 'LineWidth', linewidth,...
		'Color', '#B0C3BC' ) % 749E8E % A1BDB2
end
line([lbw hbw], [35 35], 'LineWidth', linewidth,'Color', '#014a2f')
line([CF CF], [35 75], 'LineWidth', linewidth,'Color', '#014a2f')
line([lbw lbw], [0 35.5], 'LineWidth', linewidth,'Color', '#014a2f')
line([hbw hbw], [0 35.5], 'LineWidth', linewidth,'Color', '#014a2f')

% Figure Properties
ylimits = [0 90];
xlimits = [275 4900];

ylim(ylimits)
yticks([])
xlim(xlimits)
xticks([300 fpeaks(2) 1200 fpeaks(end) 4800])
xticklabels([])
ylabel('Mag.')
set(gca,'fontsize',fontsize)
set(gca, 'XScale', 'log')
title('WB-TIN Stimulus', 'FontSize',titlesize)
xlabel('Log Freq. (Oct. w.r.t. CF)')

%% Run SFIE model for NB and WB-TIN 
clear params 

CF = 2000;
params{1}.type = 'MTFN';
params{1}.ramp_dur = 0.05;
params{1}.noise_state = 0;
params{1}.noise_band = [100, 10000];
params{1}.dur = 1;
params{1}.reptim = 1.5;
params{1}.fms = [2, 600, 3];
params{1}.mdepths = [0,0,1];
params{1}.binmode = 2;
params{1}.No = 30;
params{1}.spl = 30;
params{1}.raised_sine = 1;
params{1}.onsetWin = 25;
params{1}.Fs = 100000;
params{1}.mnrep = 3;
params{1} = generate_MTF(params{1});
params{1}.num_stim = size(params{1}.stim, 1);

params{2}.CF = CF;
params{2}.version = 5;
params{2}.dur = 0.3;
params{2}.reptim = 0.6;
params{2}.ramp_dur = 0.01;
params{2}.No = 20;
params{2}.SNR = 10:10:40;
params{2}.binmode = 2; % binaural
params{2}.nrep = 30;
params{2}.choice = 2; % Frozen noise (1), warm noise (2), all warm (3)
params{2}.condition = 0; % NoSo (0) or NoSpi (1)
params{2}.type = 'TIN_NB';
params{2}.bw_choice = 1; % NB (1) or WB(2) or 100 Hz(3)
params{2}.Fs = 100000;
params{2}.mnrep = 3;
params{2} = generate_TIN(params{2});
params{2}.num_stim = size(params{2}.stim, 1);

params{3}.type = 'WBTIN';
params{3}.version = 5;
params{3}.fpeak_mid = CF;
params{3}.ramp_dur = 0.01;
params{3}.dur = 0.3;
params{3}.reptim = 0.6;
params{3}.CF = CF;
params{3}.stp_otc = 5;
params{3}.bandwidth = 4;
params{3}.No = 20;
params{3}.SNR = [20, 30, 40];
params{3}.range = 3;
params{3}.nrep = 20;
params{3}.binmode = 2;
params{3}.onsetWin = 25;
params{3}.Fs = 100000;
params{3}.mnrep = 3;
params{3}.physio = 0;
params{3} = generate_WBTIN(params{3});
params{3}.num_stim = size(params{3}.stim, 1);

% Model parameters
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 1; % how many times to run the AN model
model_params.fiberType = 3; % AN fiber type. (1=low SR, 2=medium SR, 3=high SR)
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw 
model_params.noiseType = 0; % 0 = fixed fGn, 1 = variable fGn) 
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response
model_params.BMF = 100;
model_params.CF_range = CF;
model_params.CFs = CF;
model_params.type = 'SFIE';

% Run AN Model
num_stim = size(params, 2);
AN = cell(1, num_stim);
for ist = 1:num_stim
	timerVal = tic;
	AN{ist} = modelAN(params{ist}, model_params);
	disp(['AN model took ' num2str(toc(timerVal)/60) ' minutes'])
end

% Run SFIE model 
SFIE = cell(num_stim, 1);
for ist = 1:num_stim
	SFIE{ist} = wrapperIC(AN{ist}.an_sout, params{ist}, model_params);
end

%% Plot BE 

% MTF 
h(3) = subplot(3, 3, 3);
plotModelMTF(params{1}, SFIE{1}.average_ic_sout_BE);
yticklabels([])
set(gca, 'fontsize', fontsize)
ylim([0 40])
title({'Band-Enhanced', 'MTF Prediction'})
xlabel('Mod. Freq. (Hz)')
ylim([19 35])
yticklabels([])
ylabel('Rate')

% NB
h(4) = subplot(3, 3, 4);
plotModelTIN(params{2}, SFIE{2}.average_ic_sout_BE);
set(gca, 'fontsize', fontsize)
ylabel('Rate')
yticklabels([])
ylim([13 27])
xlabel('SNR')
xticklabels([])

% WB
colors = {'#82BB95', '#3F985C', '#034E1C'};
h(5) = subplot(3, 3, 5);
plotModelWBTIN(params{3}, SFIE{3}.average_ic_sout_BE, CF, colors);
set(gca, 'fontsize', fontsize)
ylabel('Rate')
ylim([22 28])
yticklabels([])
xlabel('Tone Freq. (kHz)')
xticklabels([])
yticklabels([])

% Plot BS 
% MTF 
h(6) = subplot(3, 3, 6);
plotModelMTF(params{1}, SFIE{1}.average_ic_sout_BS);
set(gca, 'fontsize', fontsize)
ylim([21 30])
title({'Band-Suppressed', 'MTF Prediction'})
xlabel('Mod. Freq. (Hz)')
xticklabels([])
yticklabels([])

% NB
h(7) = subplot(3, 3, 7);
plotModelTIN(params{2}, SFIE{2}.average_ic_sout_BS);
set(gca, 'fontsize', fontsize)
ylim([25 58])
xlabel('SNR')
xticklabels([])
yticklabels([])

% WB
colors = {'#82BB95', '#3F985C', '#034E1C'};
h(8) = subplot(3, 3, 8);
plotModelWBTIN(params{3}, SFIE{3}.average_ic_sout_BS, CF, colors);
set(gca, 'fontsize', fontsize)
ylim([24 32])
yticklabels([])
xlabel('Tone Freq. (kHz)')
xticklabels([])

%% Arrange 

left = [0.07 0.49 0.775];
bottom = [0.07 0.42 0.75];
width = 0.21;
height = 0.22;

set(h(1), 'position', [left(1) bottom(2) 0.32 height])
set(h(2), 'position', [left(1) bottom(1) 0.32 height])
set(h(3), 'position', [left(2) bottom(3) width 0.11])
set(h(4), 'position', [left(2) bottom(2) width height])
set(h(5), 'position', [left(2) bottom(1) width height])
set(h(6), 'position', [left(3) bottom(3) width 0.11])
set(h(7), 'position', [left(3) bottom(2) width height])
set(h(8), 'position', [left(3) bottom(1) width height])

% Create arrow
annotation('arrow',[0.249 0.349],[0.279 0.279], 'LineWidth',1, ...
	'HeadLength',5, 'HeadWidth',4);
annotation('textbox',[0.344,0.25,0.13,0.062],'String',{'Tone Freq.'},...
	'FontSize',fontsize,'EdgeColor','none');
annotation('arrow',[0.223 0.115],[0.279 0.279], 'LineWidth',1, ...
	'HeadLength',5, 'HeadWidth',4);
annotation('textbox',[0.255 0.553 0.0739 0.0623],'String',{'SNR'},...
	'FontSize',fontsize,'EdgeColor','none');
annotation('arrow',[0.252 0.252],[0.53 0.617], 'LineWidth',1, ...
	'HeadLength',5, 'HeadWidth',4);

labels = {'A', 'B', 'C', 'D'};
locs = [0.01 0.72; 0.01 0.36; 0.42 0.98; 0.71 0.98];
for ilabel = 1:4
	annotation('textbox',[locs(ilabel,:) 0.071 0.044],...
		'String',labels{1},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% Export

if save_fig == 1
	saveFigures('Fig1_hypothesis')
end
end
