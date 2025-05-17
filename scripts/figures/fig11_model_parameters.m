function fig11_model_parameters(save_fig)
% model_intermediate plots Figure 11.
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
% Created: 2025-01-21; Last revision: 2025-03-07
%
% -------------------------------------------------------------------------

%% Set up figure
[~, datapath, ~, ppi] = get_paths();

figure('Position',[50,50,3.3*ppi,3.2*ppi]);
tokensize = [10,8];
titlesize = 8;
fontsize = 7;
linewidth = 1;
labelsize = 13;

%% How strength parameter affects responses

load(fullfile(datapath, 'ModelParams1.mat'),...
	"model", "params", "CF")
h(1) = subplot(3, 2, 1);
plotThreeMTF(model, params, 0)
ylim([0 26])
xticklabels([])
set(gca, 'fontsize', fontsize)
title('MTF', 'FontSize',titlesize)

h(2) = subplot(3, 2, 2);
plotThreeWB(model, params, CF, 2)
set(gca, 'fontsize', fontsize)
xticklabels([])
leg(1) = legend(h(2), {'', '0.1', '', '0.3', '', '0.5'}, 'Location','northeastoutside');
leg(1).Title.String = {'Inhibitory', 'Strength'};
leg(1).ItemTokenSize = tokensize;
title('WB-TIN', 'FontSize',titlesize)


%% How CF range parameter affects responses
load(fullfile(datapath, 'ModelParams2.mat'),...
	"model", "params", "CF")

h(3) = subplot(3, 2, 3);
plotThreeMTF(model, params, 1)
ylim([0 18])
set(gca, 'fontsize', fontsize)
xtickangle(0)

colors2 = {'#fd8d3c','#e6550d','#a63603'};
for iWB = 2
	h(4) = subplot(3, 2, 4);
	for iparam = 1:3
		lateral_model = model{iparam}{2+iWB};

		% Analyze model
		[SNRs,~,si] = unique([params{2+iWB}.mlist.SNR].');
		num_SNRs = length(SNRs);
		[fpeaks,~,fi] = unique([params{2+iWB}.mlist.fpeak].');
		num_fpeaks = length(fpeaks);
		rate_size = [num_fpeaks,num_SNRs];
		[avIC,~,~,~] = accumstats({fi,si},lateral_model.avIC, rate_size);

		% Plot model
		hold on
		plot(fpeaks/1000, avIC(:,2), 'LineWidth',linewidth, 'Color',colors2{iparam})
		xlim([params{2+iWB}.freq_lo params{2+iWB}.freq_hi]./1000)
		set(gca, 'XScale', 'log');
		grid on
		ylim([0 32])
		set(gca, 'fontsize', fontsize)
	end
	xline(CF/1000, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
	if iWB == 1
		ylabel('Rate (sp/s)')
	end
	if iWB == 2
		xlabel('Tone Freq. (kHz)')
		set(gca, 'fontsize', fontsize)
	end
end
leg(2) = legend(h(4), {'0.25', '0.75', '1.25'}, 'Location','northeastoutside');
leg(2).Title.String = {'CF Range', '(oct)'};
leg(2).ItemTokenSize = tokensize;


%% How chosen CF affects responses
load(fullfile(datapath, 'ModelParams4.mat'),...
	"model", "CF","params1", "params2", "params3")

h(5) = subplot(3, 2, 5);
colors1 = {'#969696','#636363', '#252525'};
colors2 = {'#fd8d3c','#e6550d','#a63603'};
for iparam = 1:3
	if iparam == 1
		params = params1;
	elseif iparam == 2
		params = params2;
	else
		params = params3;
	end
	lateral_model = model{iparam}{1};
	[~, avIC, ~] = plotMTF(params{1}, lateral_model.avIC, 0);
	yline(avIC(1), 'Color',colors1{iparam}, 'LineWidth',linewidth)
	hold on
	plot(params{1}.all_fms, avIC, 'LineWidth',linewidth, 'Color', colors2{iparam})
	ylim([0,20])
	xlim([params{1}.all_fms(2) params{1}.all_fms(end)])
	set(gca,'xtick',[1.2,2, 5,  20, 50, 200],'xticklabel',...
		{'Unmod','2','5','20', '50','200'},'xscale','log')
	xtickangle(0)
	set(gca, 'XScale', 'log');
	ylabel('Rate (sp/s)','fontsize',10)
	xlabel('Mod. Freq. (Hz)')
	grid on
	set(gca, 'fontsize', fontsize)
end

for iWB = 2
	h(6) = subplot(3, 2, 6);
	for iparam = 1:3
		lateral_model = model{iparam}{2+iWB};
		if iparam == 1
			params = params1;
		elseif iparam == 2
			params = params2;
		else
			params = params3;
		end

		% Analyze model
		[SNRs,~,si] = unique([params{2+iWB}.mlist.SNR].');
		num_SNRs = length(SNRs);
		[fpeaks,~,fi] = unique([params{2+iWB}.mlist.fpeak].');
		num_fpeaks = length(fpeaks);
		rate_size = [num_fpeaks,num_SNRs];
		[avIC,~,~,~] = accumstats({fi,si},lateral_model.avIC, rate_size);

		% Plot model
		fpeaks_re_CF = log2(fpeaks./CF(iparam));
		yline(avIC(1,1), 'k', 'LineWidth',linewidth, 'Color', colors1{iparam})
		hold on
		plot(fpeaks_re_CF, avIC(:,2), 'LineWidth',linewidth, 'Color', colors2{iparam})
		yticks([0 10 20 30])
		if iWB == 1
			ylabel('Rate (sp/s)')
		elseif iWB == 2
			xlabel('Tone Freq. w.r.t. CF (Oct)')
		end
	end
	xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
	ylim([0,32])
	xlim([-1.5 1.5])
	grid on
	set(gca, 'fontsize', fontsize)
end
leg(3) = legend(h(6), {'', '1000', '', '3000', '', '5000'}, 'Location','northeastoutside');
leg(3).Title.String = 'CF (Hz)';
leg(3).ItemTokenSize = tokensize;

%% Arrange Figure

height = 0.22;
width = [0.3 0.29 0.3 0.29 0.3 0.29];
left = [0.12 0.5];
bottom = [0.72, 0.44, 0.09];

row = reshape(repmat(1:3, 2, 1), 6, 1);
col = repmat(1:2, 1, 3);
for ii = 1:6
	irow = row(ii);
	icol = col(ii);
	set(h(ii), 'position', [left(icol) bottom(irow) width(ii) height])
end

leg(1).Position = [0.876, 0.79, 0.0757, 0.11];
leg(2).Position = [0.876, 0.52, 0.0757, 0.11];
leg(3).Position = [0.876, 0.19, 0.0757, 0.11];

labels = {'A', 'B', 'C'};
labelbottom = [0.975 0.67 0.33];
for ii = 1:3
	annotation('textbox',[0.01 labelbottom(ii) 0.071 0.044],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% Export

if save_fig == 1
	saveFigure('Fig11_model_parameters')
end
end