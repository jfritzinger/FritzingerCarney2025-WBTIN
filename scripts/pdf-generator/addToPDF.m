function [img, images] = addToPDF(cluster, params, rabbit, images, fig)
import mlreportgen.dom.*

% Set figure size
if iscell(params)
	param = params{1}; % If multiple, plots only the first instance
else
	param = params;
end

type = param.type;
switch param.type
	case 'char_spl'
		values = [2.7 2];
		param.binmode = 2;
	case 'typMTFN'
		values = [2.5 2];
	case 'type=RM'
		type = 'RM';
		values = [];
	case 'SCHR'
		values = [7.5 3];
	case 'STRF'
		values = [6 2.5];
		type = param.plot_type;
	case 'SPEC_slide'
		type = param.plot_type;
		switch param.SPEC_slide_type
			case 'WB_noise'
				type = param.plot_type;
				if strcmp(param.plot_type, 'WBTIN') ||...
						strcmp(param.plot_type, 'WBTIN_Onset')||...
						strcmp(param.plot_type, 'WBTINRM_Comparison')||...
						strcmp(param.plot_type, 'WBTINneuro')
					values = [7.5 2.5*ceil(param.num_DSIDs/3)];
				elseif strcmp(param.plot_type, 'WBTINRasters')
					values = [8.5 param.num_plots*0.17];
				elseif strcmp(param.plot_type, 'WBTINPSTH')
					values = [8.5 param.num_plots*0.17];
				elseif strcmp(param.plot_type, 'WBTIN_onCF')
					values = [4 2.5];
				end
			case 'NB_noise'
				type = param.plot_type;
				if strcmp(param.plot_type, 'NBTIN')||...
						strcmp(param.plot_type,'NBTIN_subnoise')
					if param.version < 5
						values = [7.5 3];
					else
						values = [7.5 2.5*ceil(param.num_DSIDs/3)];
					end
				elseif strcmp(param.plot_type, 'NBTIN_onCF')
					values = [4 2.5];
				end
		end
end

if ~strcmp(param.type, 'type=RM')
	changeFigureSize(values)
end


% Add the plot to the document
name = sprintf('%s_%s_TT%d_N%d_%d_%d.svg', type, rabbit, cluster.tetrode,...
	cluster.neuron, param.dsid, param.binmode);
tempfilepath = name;
print(fig, tempfilepath, '-dsvg');
img = Image(tempfilepath);
delete(fig) %delete plot figure window
images = [images {img}];

end
