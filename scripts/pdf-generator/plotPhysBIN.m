function [fig, data] = plotPhysBIN(cluster, params, stim)
% plotPhysBIN plots the average rate response to binaural noise.
%   Average rate response for binaural, contralateral, and ipsilateral
%   noise that is presented as DSID 1 in all physiology recordings.
%
% Inputs:
%    cluster - struct, output of post_process
%    params - params from post_process that have the type 'char_spl'
%	 stim - struct, output of post_process
%
% Outputs:
%    fig - Average rate plot
%    data - struct that contains processed average rate data. 
%		 data.rate - average rate response (5x3, num_spls x num_binmodes)
%		 data.rate_std - standard deviation of rates (5x3, num_spls x num_binmodes)
%		 data.spls - all spls used (5)
%		 data.num_binmodes - 3, binaural, contralateral, and ipsilateral
%
% Other m-files required: accumstats, checkIfSPLShift, isDriven
%
% Author: J. Fritzinger, adapted from cell_characterization
% Created: 2022-06-18; Last revision: 2024-09-08
%
% -------------------------------------------------------------------------

% Process cluster 
if iscell(params)
	param = params{1}; % If multiple, plots only the first instance
else
	param = params;
end
ds = param.dsid;
dsi = 1;
this_ds = stim(ds).dsid == ds(dsi);

% Spike rate vs. SPL
[spls,~,si] = unique([param.list.spl].');
[binmodes,~,bmi] = unique([param.list.binmode].');
num_spls = length(spls);
num_binmodes = length(binmodes);
rate_size = [num_spls,num_binmodes];
dur = param.dur/1000; % stimulus duration in seconds
num_spikes = cluster.num_spikes_peri(this_ds);
[rate,rate_std] = accumstats({si,bmi},num_spikes/dur,rate_size); 

% Check if this session has an spl shift 
spl_shift = checkIfSPLShift(param.animal, param.session);
spls = spls + spl_shift;
param.spls = spls;

% Plot
fig = figure;
box on
hold on
h = errorbar(repmat(spls,1,num_binmodes),rate,rate_std);
hold off
set(h,{'Color'},{'r';[0 0.8 0];'k'})
hLegend = legend(h,{'Ipsilateral','Contralateral','Binaural'},...
	'Location','Northeast', 'FontSize',3);
hLegend.ItemTokenSize = [6,6];
cur_ax = gca; % save for later plotting
set(cur_ax,'FontSize',7)

if min(spls) < 0 % for the new setup
	isspont = 1;
	xlimits = [spls(2),spls(end)] + [-2 2];
elseif max(spls) > 30 % old Characterize miswrote param of spont spl (-100 dB) as the default spl (65 dB), but correct with the stimuli.
	isspont = 1;
	xlimits = [spls(1),spls(end-1)] + [-2 2]; %  Before 2015, the spont rate was not visible due to this xlim
else
    xlimits = [spls(1),spls(end)]+[-2 2];
	isspont = 0;
end

xlim(xlimits)
ylimits = ylim;
span = diff(ylimits);
new_ylim = ylimits + [-span 0]/50;
ylim(new_ylim)
xlabel('Noise Spectrum Level (dB SPL)')
ylabel('Spike Rate (sp/s)')
title(['Binaural Response, Onset Win = ' num2str(param.onsetWin) 'ms'])

% Add line at spontaneous rate. Referenced to MTFN analysis to add a line.
% Will have errors with the early new setup session because the spont rate
% was not obtained yet.
if isspont == 1
	hold on;
	spont = sum(num_spikes(1:3))/(3*dur); % 3 is hard-coded because the # of rep is 3.
	line(xlimits,[1 1]*spont,'Color',[0.5 0.5 0.5])
	% Move spontaneous rate line behind other lines.
	kids = get(gca,'Children');
	set(cur_ax,'Children',kids([2:end,1]))
	hold off;
end

driven = isDriven(cluster, params, stim);
hold on
if driven == 1
     dim = [0.15 0.6 0.3 0.3];
    annotation('textbox',dim,'String','Driven','FitBoxToText','on');
else
     dim = [0.15 0.6 0.3 0.3];
    annotation('textbox',dim,'String','NOT DRIVEN','FitBoxToText','on');
end

% Create struct that contains processed data 
data.rate = rate;
data.rate_std = rate_std;
data.spls = spls;
data.num_binmodes = num_binmodes;

end