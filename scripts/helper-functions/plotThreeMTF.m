function plotThreeMTF(model, params, labels)
%PLOTTHREEMTF Plots three Modulation Transfer Functions (MTFs) for model comparison
%   plotThreeMTF(MODEL, PARAMS, LABELS) visualizes three MTF curves with 
%   standardized formatting for model evaluation.
%
%   Inputs:
%       model  - Cell array {3×1} containing model structures with fields:
%                .avIC - Average inferior colliculus response rates
%       params - Cell array {1×1} containing parameters structure with:
%                .all_fms - Modulation frequencies (Hz)
%                .mnrep   - Number of repetitions (for potential error calculations)
%       labels - Binary flag (0/1) to control x-axis labeling

linewidth = 1;
colors1 = {'#969696','#636363', '#252525'};
colors2 = {'#fd8d3c','#e6550d','#a63603'};
for iparam = 1:3
	lateral_model = model{iparam}{1};
	[~, avIC, ~] = plotMTF(params{1}, lateral_model.avIC, 0);
	yline(avIC(1), 'Color',colors1{iparam}, 'LineWidth',linewidth)
	hold on
	%errorbar(params{1}.all_fms, avIC, stdIC./sqrt(params{1}.mnrep), 'LineWidth',linewidth)
	plot(params{1}.all_fms, avIC, 'LineWidth',linewidth, 'Color',colors2{iparam})
	xlim([params{1}.all_fms(2) params{1}.all_fms(end)])
	set(gca,'xtick',[1.2,2, 5,  20, 50, 200],'xticklabel',...
		{'Unmod','2','5','20', '50','200'},'xscale','log')
	set(gca, 'XScale', 'log');
	ylabel('Rate (sp/s)','fontsize',10)
	grid on
	if labels == 1
		xlabel('Mod. Freq. (Hz)')
	end
end
end