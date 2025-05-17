%% Fig0_PutativeNeuronPDF
% J. Fritzinger, updated 10/27/23
%
% Reads in spreadsheet and creates PDF for all sessions (sorted by CF?)
clear 

import mlreportgen.dom.*
import mlreportgen.report.*
timerVal = tic;


%% Load in spreadsheet

% Set paths & load in spreadsheet
[base, datapath, savepath] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
filename = 'PutativeNeuronPDF_generated';

%% Set up PDF report if needed

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
report_name = fullfile(savepath, sprintf('%s.pdf', filename));
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.01in';
pm.PageMargins.Header = '0.01in';
pm.PageMargins.Bottom = '0.01in';
pm.PageMargins.Footer = '0.01in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

% Sort the CFs by their indices on the spreadsheet
CFs = sessions.CF;
[~, order] = sort(CFs);

%% Create PDF

num_neurons = size(sessions,1);
for c_ind = 1 %1:num_neurons % Loops through each cluster

	% Label the session
	putative_neuron = sessions.Putative_Units{order(c_ind)};
	MTF_shape = sessions.MTF{order(c_ind)};
	CF = sessions.CF(order(c_ind));
    % putative_neuron = sessions.Putative_Units{c_ind};
	% MTF_shape = sessions.MTF{c_ind};
	% CF = sessions.CF(c_ind);

	% Load in data
	load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');

	% Create heading
	h = Heading(2, sprintf('%s, CF = %0.0fHz, %s', ...
		putative_neuron,  CF, MTF_shape));
	b = Border();
	b.BottomStyle = 'single';
	b.BottomColor = 'LightGray';
	b.BottomWidth = '1pt';
	h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	append(rpt,h);
	fprintf('Creating plots... %s, CF = %0.0fHz\n', putative_neuron, CF);

	% Binaural Response
	% params_BIN = data(1,1); % Binaural noise response
	% yes_BIN = ~cellfun(@isempty,params_BIN);
	% if any(yes_BIN)
	% 	[fig, data_BIN] = plotPhysBIN(params_BIN{1}.cluster, params_BIN, params_BIN{1}.stim);
	% 	[BINplt, images] = addToPDF(params_BIN{1}.cluster, params_BIN, putative_neuron, images, fig);
	% end

	for ibin = 2

		% Get binaural data for each stimulus
		params_RM = data(2,ibin); % Gets response map
		params_MTF = data(3,ibin); % Gets MTFN
		params_STRF = data(4,ibin);
		params_WB = data(6:8,ibin); % Gets WB-TIN
		params_NB = data(9:11,ibin); % Gets NB-TIN
		params_ITD = data(12:end,1);


		% Response Map
		yes_RM = ~cellfun(@isempty,params_RM);
		if any(yes_RM)
			[fig, data_RM] = plotPhysRM(params_RM{1}.cluster, params_RM,params_RM{1}.stims, CF);
			[rmplt, images] = addToPDF(params_RM{1}.cluster, params_RM, putative_neuron, images, fig);
			rmplt.Height = '2.6in';
			rmplt.Width = '5.3in';
		end

		% MTF
		yes_MTF = ~cellfun(@isempty,params_MTF);
		if any(yes_MTF)
			[~,~,~, fig, data_MTF] = plotPhysMTF_TTest(params_MTF{1}.cluster, params_MTF,params_MTF{1}.stims);
			[mtfplt, images] = addToPDF(params_MTF{1}.cluster, params_MTF, putative_neuron, images, fig);
		end

		if any(yes_RM) && any(yes_MTF)
			rm_mtf_table = Table({rmplt,mtfplt});
			rm_mtf_table.Style = {Width('100%'), ResizeToFitContents(false)};
			rm_mtf_table.BorderColor = 'White';
			append(rpt, rm_mtf_table);
		elseif any(yes_RM)
			append(rpt, rmplt);
		elseif any(yes_MTF)
			append(rpt, mtfplt);
		end

		% WB-TIN
		yes_WB = ~cellfun(@isempty,params_WB);
		num_DSIDs = sum(yes_WB);
		params_WB = params_WB(yes_WB);
		data_WBTIN = cell(num_DSIDs, 1);
		if any(yes_WB)
			[params_WB, fig, data_WBTIN] = plotPhysWBTIN(params_WB{1}.cluster, params_WB, [], CF, data_RM, data_WBTIN, 'No');
			[wbplt, images] = addToPDF(params_WB{1}.cluster, params_WB, putative_neuron, images, fig);
			append(rpt, wbplt);

			% On-CF
			[params_WB, fig, data_WBTIN] = plotPhysWBTIN_onCF(params_WB{1}.cluster, params_WB, [], CF, data_WBTIN);
			[CFwbplt, images] = addToPDF(params_WB{1}.cluster, params_WB, putative_neuron, images, fig);

		end

		% NB-TIN
		%yes_NB = [];
		yes_NB = ~cellfun(@isempty,params_NB);
		num_DSIDs = sum(yes_NB);
		params_NB = params_NB(yes_NB);
		data_NBTIN = cell(num_DSIDs, 1);
		if any(yes_NB)
			[params_NB, fig, data_NBTIN] = plotPhysNBTIN(params_NB{1}.cluster, params_NB, [], CF, 'Normal');
			[nbplt, images] = addToPDF(params_NB{1}.cluster, params_NB, putative_neuron, images, fig);
			append(rpt, nbplt);

			% On-CF
			[params_NB, fig, data_NBTIN] = plotPhysNBTIN_onCF(params_NB{1}.cluster, params_NB, [], CF, data_NBTIN);
			[CFnbplt, images] = addToPDF(params_NB{1}.cluster, params_NB, putative_neuron, images, fig);
		end

		% Plots WBTIN and NBTIN on-CF plots together if possible
		if any(yes_NB) && any(yes_WB)
			onCF_table = Table({CFnbplt,CFwbplt});
			onCF_table.Style = {Width('100%'), ResizeToFitContents(false)};
			onCF_table.BorderColor = 'White';
			append(rpt, onCF_table);
		elseif any(yes_NB)
			append(rpt, CFnbplt);
		elseif any(yes_WB)
			append(rpt, CFwbplt);
		end

	end
end


% Close the report and delete all temporary files
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)

% Time how long the PDF creation took
elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])