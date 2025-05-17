function [params, fig, data] = plotPhysWBTIN_onCF(cluster, params, stims, CF, data)
% Plot the on-CF response to the sliding WB-TIN

num_DSIDs = length(params);
binmode_strs = {'Ipsi','Contra','Bin'};
fig = figure;
for ind = 1:num_DSIDs

    % Parameters
    param = params{ind};
    ds = param.dsid;
	if isempty(stims)
		this_ds = param.stims.dsid == ds;
		cluster = param.cluster;
	else
		this_ds = stims(ds).dsid == ds;
	end
    binmode = binmode_strs{param.binmode + 1};

    % Check if this session has an spl shift
    spl_shift = checkIfSPLShift(param.animal, param.session);
    No = param.No + spl_shift;

    % Analysis
    if param.version <= 2 % Error in stimulus
        h = subplot(1, 1, 1);
        pos = get(h, 'position');
        dim = pos.*[1 1 0.5 0.5];
        str =  ['This presentation was version ' num2str(param.version)];
        annotation('textbox',dim, 'String',str, 'FitBoxToText','on', 'verticalalignment', 'bottom')
        set(gca,'FontSize',7)
		params{1}.plot_type = 'WBTIN_onCF';
		params{1}.num_DSIDs = num_DSIDs;
		data{ind}.rate_onCF = [0 0 0];
		data{ind}.rate_std_onCF = [0 0 0];
		data{ind}.p_onCF = [];
		data{ind}.CF = [];
	elseif isnan(CF)
		h = subplot(1, 1, 1);
		pos = get(h, 'position');
		dim = pos.*[1 1 0.5 0.5];
		str =  'No CF Recorded';
		annotation('textbox',dim, 'String',str, 'FitBoxToText','on', 'verticalalignment', 'bottom')
		set(gca,'FontSize',7)
		params{1}.plot_type = 'WBTIN_onCF';
		params{1}.num_DSIDs = num_DSIDs;
		box on
		title('On-CF WBTIN')
		data{ind}.rate_onCF = [0 0 0];
		data{ind}.rate_std_onCF = [0 0 0];
		data{ind}.p_onCF = [];
		data{ind}.CF = [];
	else
		[fpeaks,~,fi] = unique([param.list.fpeak].');
		num_fpeaks = length(fpeaks);
        [SNR,~,si] = unique(double([param.list.SNR]).');
        SNR = flip(SNR);
        num_SNRs = length(SNR);

        if param.dur > 1
            dur = param.dur/1000; % stimulus duration in seconds.
        else
            dur = param.dur;
        end

        rate_size = [num_fpeaks,num_SNRs];
        spike_rates = cluster.num_spikes_delayed(this_ds)/...
            (dur - param.onsetWin/1000);


        if length(fi) == length(spike_rates)
            [rate,rate_std,~,~] = accumstats({fi,si},spike_rates, rate_size);

            % Find rate nearest to CF
            SNRs(ind,1) = categorical(-inf);
            rate_onCF(ind,1) = mean(rate(1,:));
            rate_std_onCF(ind,1) = mean(rate_std(1,:));

            SNRs(ind,2:num_SNRs+1) = categorical(param.SNR);
            [~,closestIndex] = min(abs(fpeaks(:, 1)-CF));
            rate_onCF(ind, 2:num_SNRs+1) = rate(closestIndex,:);
            rate_std_onCF(ind,2:num_SNRs+1) = rate_std(closestIndex,:);
			if sum(rate_onCF(ind,:))>0
				ttest_r = ttest_meanandstd(rate_onCF(ind,1), rate_std_onCF(ind,1), ...
					rate_onCF(ind,num_SNRs+1), rate_std_onCF(ind,num_SNRs+1), param.nrep);
				if ttest_r<0.05
					label{ind} = sprintf('DS%d, No=%d, p = %.3f*', ds, No, ttest_r);
				else
					label{ind} = sprintf('DS%d, No=%d, p = %.3f', ds, No, ttest_r);
				end
			else
				label{ind} = 'Reanalyze, error';
				ttest_r = [];
			end
			data{ind}.rate_onCF = rate_onCF(ind, :);
			data{ind}.rate_std_onCF = rate_std_onCF(ind, :);
			data{ind}.p_onCF = ttest_r;
			data{ind}.CF = CF;
		else
			data{ind}.rate_onCF = [0 0 0];
			data{ind}.p_onCF = [];
			data{ind}.CF = [];
			data{ind}.rate_std_onCF = [0 0 0];
			label{ind} = 'Reanalyze, error';
			ttest_r = [];
        end
    end
    clear leg
end

% Plot on-CF results
if param.version > 2 && ~isnan(CF)
	errorbar(SNRs', rate_onCF', (rate_std_onCF/sqrt(param.nrep))', 'Linewidth', 2);
	hLegend = legend(label, 'location', 'best', 'FontSize',6);
	hLegend.ItemTokenSize = [6,6];
	xlabel('SNR (dB SPL)')
	ylabel('Avg. Rate (sp/s')
	grid on
	box on
	set(gca,'FontSize',7)
	title(sprintf('On-CF WBTIN, CF=%dHz, Stim=%dHz', round(CF), round(fpeaks(closestIndex))));
	axis_set = ylim;
	axis_set(1) = 0;
	ylim(axis_set);
end

% Add type of plot
params{1}.plot_type = 'WBTIN_onCF';
params{1}.num_DSIDs = num_DSIDs;
end

function tprob = ttest_meanandstd(A, A_std, B, B_std, N)

v = 2*N-2;
tval = (A-B) / sqrt((A_std^2+B_std^2)/N);       % Calculate T-Statistic
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
%tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
%tprob = 1-[tdist2T(tval,v)  tdist1T(tval,v)];
tprob = 1-tdist2T(tval,v);

end