function avIC = plotModelMTF(params, IC_response)
% Visualizes MTF model responses
%   AVIC = PLOTMODELMTF(PARAMS, IC_RESPONSE) generates a MTF plot from model
%   responses and returns average IC responses.
%
%   Inputs:
%       params - Struct containing modulation parameters:
%           .all_fms : Modulation frequencies (Hz)
%       IC_response - Model response data array
%
%   Output:
%       avIC - Average response vector

linewidth = 1;

spont_color = [0.4 0.4 0.4];
[~, avIC, ~] = plotMTF(params, IC_response, 0);
yline(avIC(1), 'Color',spont_color, 'LineWidth',linewidth)
hold on
plot(params.all_fms, smooth(avIC), 'LineWidth',linewidth,'Color','k')
xlim([params.all_fms(2) params.all_fms(end)])
set(gca,'xtick',[1.2,2, 5,  20, 50, 200, 500],'xticklabel',...
	{'Unmod','2','5','20', '50','200','500'},'xscale','log')
set(gca, 'XScale', 'log');
xticklabels([])

end
