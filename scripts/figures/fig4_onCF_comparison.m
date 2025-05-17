function fig4_onCF_comparison(save_fig)
% onCF_comparison plots Figure 3. 
%	This function plots scatter plots of on-CF WB vs NB TIN rates and
%	separates the data by BE or BS MTF types
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN
%	 .mat files required: -
%	 spreadsheets required: Data_Table.xlsx
%
% Author: J. Fritzinger
% Created: 2025-01-20; Last revision: 2025-03-06
%
% -------------------------------------------------------------------------


%% Load in data

[~, datapath, ~, ppi] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure

figure('Name', 'On-CF WB NB Comparison', 'Position', [50,50,6.929*ppi,2.4*ppi]); % left bottom width height

subplot_numbers = [1, 4, 7; 10, 13, 16];
MTF_colors = {'b', 'r', '#7F7FFD', '#FF7E7E', '#77AC30'};
legsize = 6;
titlesize = 8;
fontsize = 7;
labelsize = 13;

%% Find sessions of interest

WB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_3);
WB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_23);
WB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_43);
WB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_3);
WB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_23);
WB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_43);

NB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_con_3);
NB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_con_23);
NB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_con_43);
NB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_3);
NB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_23);
NB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_NB_noise_43);

contra = WB_noise_con & NB_noise_con;
bin = WB_noise & NB_noise;

%% For each session:

binmodes = {'Contra', 'Binaural'};
Nos = {'3', '23', '43'};
for ibin = 2 
	for iNo = 1:3 % For each No

		iBE = 1; iBS = 1; iH = 1; iF = 1; iHBE = 1; iHBS = 1;
		MTF_BE = zeros(1, 4); MTF_BS = zeros(1, 4); MTF_H = []; MTF_F = [];
		MTF_HBE = []; MTF_HBS = [];
		if ibin == 1
			index = find(contra(:,iNo));
		else
			index = find(bin(:,iNo));
		end

		for isesh = 1:length(index)
			s_ind = index(isesh);

			% Load in session
			putative_neuron = sessions.Putative_Units{s_ind};
			CF = sessions.CF(s_ind);
			load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');

			% Get data for each stimulus
			params_MTF = data{3,ibin}; % Gets binaural MTFN
			params_WB = data(5+iNo,ibin); % Gets binaural WB-TIN stimuli
			params_NB = data(8+iNo,ibin);

			if isempty(params_MTF) && ibin == 1
				params_MTF = data{3,2};
			elseif isempty(params_MTF) && ibin == 2
				params_MTF = data{3,1};
			end

			% RM and MTF analysis
			data_MTF = analyzeMTF(params_MTF);

			% Get parameters and data for each session
			data_WB = analyzeWBTIN(params_WB, CF); 
			data_WB = data_WB{1};
			data_NB = analyzeNBTIN(params_NB, CF);
			data_NB = data_NB{1};
	
			% Get change in rate for WB and NB (raw rates)
			change_nb = data_NB.rate_onCF(3) - data_NB.rate_onCF(1);
			change_wb = data_WB.rate_onCF(3) - data_WB.rate_onCF(1);
			std_nb = (sqrt(data_NB.rate_std_onCF(3)^2 + data_NB.rate_std_onCF(1)^2))/sqrt(params_NB{1}.nrep); % Propogation of uncertainty
			std_wb = (sqrt(data_WB.rate_std_onCF(3)^2 + data_WB.rate_std_onCF(1)^2))/sqrt(params_WB{1}.nrep); % Propogation of uncertainty

			all(isesh, :) = [change_nb change_wb std_nb std_wb];
			if strcmp(data_MTF.MTF_shape, 'BE')
				MTF_BE(iBE,:) = [change_nb change_wb std_nb std_wb];
				iBE = iBE + 1;
			elseif strcmp(data_MTF.MTF_shape, 'BS')
				MTF_BS(iBS,:) = [change_nb change_wb std_nb std_wb];
				iBS = iBS + 1;
			elseif contains(data_MTF.MTF_shape, 'H')
				MTF_H(iH,:) = [change_nb change_wb std_nb std_wb];
				iH = iH + 1;
				if strcmp(data_MTF.at_100, 'BE')
					MTF_HBE(iHBE,:) = [change_nb change_wb std_nb std_wb];
					iHBE = iHBE + 1;
				else
					MTF_HBS(iHBS,:) = [change_nb change_wb std_nb std_wb];
					iHBS = iHBS + 1;
				end
			elseif strcmp(data_MTF.MTF_shape, 'F')
				MTF_F(iF,:) = [change_nb change_wb std_nb std_wb];
				iF = iF + 1;
			end
		end

		h(subplot_numbers(1, iNo)) = subplot(6, 3, subplot_numbers(1, iNo)); % [1 4 7 10 13 16]
		h(subplot_numbers(1, iNo)+1) = subplot(6, 3, subplot_numbers(1, iNo)+1);
		h(subplot_numbers(1, iNo)+2) = subplot(6, 3, subplot_numbers(1, iNo)+2);
		for iMTF = 1:4
			if iMTF == 1
				nbwb_data = MTF_BE;
			elseif iMTF == 2
				nbwb_data = MTF_BS;
			elseif  iMTF == 3
				nbwb_data = MTF_HBE;
			elseif iMTF == 4
				nbwb_data = MTF_HBS;
			else
				nbwb_data = MTF_F;
			end
			outliers = findOutliers(nbwb_data);
			nbwb_data(outliers,:) = [];

			% Plot 1
			a = subplot_numbers(1, iNo);
			color = MTF_colors{iMTF};
			hold(h(a), "on")
			scatter(h(a), nbwb_data(:,1), nbwb_data(:,2), 2, 'filled', 'MarkerEdgeColor',color, 'MarkerFaceColor',color)
			errorbar(h(a), nbwb_data(:,1), nbwb_data(:,2), nbwb_data(:,4)./2, ...
				nbwb_data(:,4)./2, nbwb_data(:,3)./2, nbwb_data(:,3)./2,...
				'o','color', color, 'MarkerEdgeColor',color, 'MarkerFaceColor',...
				color, 'MarkerSize',1, 'CapSize',1)		
			if iMTF == 1 || iMTF == 2
				e = calculateEllipseSTD(nbwb_data(:,1:2));
				plot(h(a), e(1,:), e(2,:), 'Color',MTF_colors{iMTF}, 'linewidth', 1);
			end

			% Create distribution plot for the X-axis (horizontal)
			b = subplot_numbers(1, iNo)+1;
			edges = -60:5:60;
			hold(h(b), "on")
			xline(h(b), 0, 'k')
			histogram(h(b), nbwb_data(:, 1),edges, 'Orientation', 'vertical', 'EdgeColor', 'k', 'FaceColor', color);
			
			% Create distribution plot for the Y-axis (vertical) below the scatter plot
			c = subplot_numbers(1, iNo)+2;
			hold(h(c), "on")
			yline(h(c), 0, 'k')
			histogram(h(c), nbwb_data(:, 2), edges,  'Orientation', 'horizontal', 'EdgeColor', 'k', 'FaceColor', color);
		end

		% Figure properties, plot a
		xlim(h(a), 60*[-1 1])
		ylim(h(a), 60*[-1 1])
		yline(h(a), 0, 'k')
		xline(h(a), 0, 'k')
		xticklabels(h(a), [])
		yticklabels(h(a), [])
		title(h(a), sprintf('N_0 = %s dB SPL  ', Nos{iNo}), 'fontsize', titlesize)
		if iNo == 1
			hLegend = legend(h(a), 'BE','','', 'BS', '','','H_B_E', '','H_B_S', 'fontsize',...
				legsize, 'location', 'southeast');
			hLegend.ItemTokenSize = [6,6];
		end

		% Figure properties, plot b
		xlabel(h(b), 'NB: r_{40}-r_{noise} (sp/s)')
		xlim(h(b), 60*[-1 1])
		xticks(h(b), [-50 -25 0 25 50])
		yticks(h(b), [0 5 10])
		yticklabels(h(b), [])
		set(h(b), 'fontsize', fontsize)

		% Figure properties, plot c
		if iNo == 1
			ylabel(h(c), 'WB:, r_{40}-r_{noise} (sp/s)')
		end
		ylim(h(c), 60*[-1 1])
		yticks(h(c), [-50 -25 0 25 50])
		xticks(h(c), [0 5 10])
		xticklabels(h(c), [])
		set(h(c), 'fontsize', fontsize)

		% Numbers of neurons per category
		outliers = findOutliers(MTF_BE);
		MTF_BE_ad = MTF_BE;
		MTF_BE_ad(outliers,:) = [];
		outliers = findOutliers(MTF_BS);
		MTF_BS_ad = MTF_BS;
		MTF_BS_ad(outliers,:) = [];

		signed_vals = [sign(MTF_BE_ad(:,1:2)); sign(MTF_BS_ad(:,1:2))];
		num_quad1 = sum(signed_vals(:,1)>0 & signed_vals(:,2)>0); 
		num_quad2 = sum(signed_vals(:,1)<0 & signed_vals(:,2)>0); 
		num_quad3 = sum(signed_vals(:,1)<0 & signed_vals(:,2)<0); 
		num_quad4 = sum(signed_vals(:,1)>0 & signed_vals(:,2)<0); 
		xL=xlim(h(a));
		yL=ylim(h(a));
		text(h(a), 0.95*xL(1),0.99*yL(2),['n=' num2str(num_quad2)],...
			'HorizontalAlignment','left','VerticalAlignment','top', 'FontSize',fontsize)
		text(h(a), 0.95*xL(2),0.99*yL(2),['n=' num2str(num_quad1)],...
			'HorizontalAlignment','right','VerticalAlignment','top', 'FontSize',fontsize)
		text(h(a), 0.95*xL(1),0.99*yL(1),['n=' num2str(num_quad3)],...
			'HorizontalAlignment','left','VerticalAlignment','bottom', 'FontSize',fontsize)
		text(h(a), 0.95*xL(2),0.99*yL(1),['n=' num2str(num_quad4)],...
			'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize',fontsize)

		% Print out percentages 
		increase = sum(all(:,1)>0);
		decrease = sum(all(:,1)<0);
		fprintf('NB: For %s dB, %s: %d/%d increase, %0.1f percent \n', ...
			Nos{iNo}, binmodes{ibin},increase, size(all, 1), increase/size(all, 1)*100);
		fprintf('NB: For %s dB, %s: %d/%d decrease, %0.1f percent \n', ...
			Nos{iNo}, binmodes{ibin},decrease, size(all, 1), decrease/size(all, 1)*100);
		increase = sum(all(:,2)>0);
		decrease = sum(all(:,2)<0);
		fprintf('WB: For %s dB, %s: %d/%d increase, %0.1f percent \n', ...
			Nos{iNo}, binmodes{ibin},increase, size(all, 1), increase/size(all, 1)*100);
		fprintf('WB: For %s dB, %s: %d/%d decrease, %0.1f percent \n \n', ...
			Nos{iNo}, binmodes{ibin},decrease, size(all, 1), decrease/size(all, 1)*100);
	end
end

%% Position plots

all_fig_positions = ...
   [0.12,0.27,0.21,0.6;...
	0.43,0.27,0.21,0.6;...
	0.74,0.27,0.21,0.6]; % left bottom width height

subplot_numbers = [1, 4, 7];
for ipos = 1:3
	fig_position = all_fig_positions(ipos,:);
	nb_position = [fig_position(1),fig_position(2)-0.12,fig_position(3),0.09];
	wb_position = [fig_position(1)-0.04,fig_position(2),0.03,fig_position(4)];
	set(h(subplot_numbers(ipos)), 'Position', fig_position)
	set(h(subplot_numbers(ipos)+1), 'Position', nb_position)
	set(h(subplot_numbers(ipos)+2), 'Position', wb_position)
end

% Create textbox
labels = {'A', 'B', 'C'};
labelleft= repmat([0.05, 0.36, 0.665], 1, 2);
labelbottom = [repmat(0.96,1, 3) repmat(0.48, 1, 3)];
for ii = 1:3
	annotation('textbox',[labelleft(ii) labelbottom(ii) 0.071 0.058],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% Export

if save_fig == 1
	saveFigure('Fig4_onCF_comparison')
end
end
