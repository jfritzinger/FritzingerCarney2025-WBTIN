function fig = compareMultiple(names, titlename, f, legsize, loc, varargin)
colors = {'b', 'r', 'g', 'm', 'c'};

for itype = 1:length(names)

	WB_array_z = squeeze(varargin{itype}(:,1, :,:));

	oct_mat = WB_array_z;
	num_samples = sum(squeeze(sum(~isnan(oct_mat),2, 'omitnan')), 1, 'omitnan');
	oct_avg_all1 = mean(oct_mat, 1, 'omitnan');
	oct_avg_all = mean(squeeze(oct_avg_all1), 1, 'omitnan');
	fprintf('%s, %s: n = %d\n', titlename, names{itype}, size(oct_avg_all1, 2))
	n(itype) = size(oct_avg_all1, 2);

	not_nan_ind = ~isnan(oct_avg_all);
	oct_avg = oct_avg_all(not_nan_ind);

	if ndims(oct_mat) == 3
		oct_std_all = squeeze(std(oct_mat,0,[1 2], 'omitnan'))./sqrt(num_samples');
		oct_std = oct_std_all(not_nan_ind)';
	elseif ndims(oct_mat) == 2
		oct_std_all = std(oct_mat,0,1, 'omitnan')./sqrt(num_samples');
		oct_std = oct_std_all(not_nan_ind);
	end
	f_notnan = f(not_nan_ind);

	hold on
	plot(f_notnan, oct_avg, 'color', 'k', 'LineWidth',1);
	patch([f_notnan flip(f_notnan)], [oct_avg-oct_std flip(oct_avg+oct_std)], ...
		colors{itype}, 'FaceAlpha',0.25, 'EdgeColor','none')

end
yline(0)
xline(0)
xlim([-1.7 1.7])
ylim([-1.5 2.5])
for ii = 1:length(names)
	leg{2*ii-1} = '';
	leg{2*ii} = sprintf('%s, n=%d', names{ii}, n(ii));
end
set(gca, 'FontSize', 12)


CF_tests = -1:0.25:1;
for ii = 1:9
	xline(CF_tests(ii), '--', 'color', [0.5 0.5 0.5])
	leg{length(names)*2+ii} = '';
end
set(gca,'XGrid','off','YGrid','on')
hleg = legend(leg, 'Location',loc, 'FontSize',legsize);
hleg.ItemTokenSize = [10,8];

end

% function fig = compareMultiple(names, titlename, f, varargin)
% colors = {'b', 'r', 'g', 'm', 'c'};
% 
% for itype = 1:length(names)
% 
% 	WB_array_z = squeeze(varargin{itype}(:,1, :,:));
% 
% 	oct_mat = WB_array_z;
% 	num_samples = sum(squeeze(sum(~isnan(oct_mat),2, 'omitnan')), 1, 'omitnan');
% 	oct_avg_all1 = mean(oct_mat, 1, 'omitnan');
% 	oct_avg_all = mean(squeeze(oct_avg_all1), 1, 'omitnan');
% 	fprintf('%s, %s: n = %d\n', titlename, names{itype}, size(oct_avg_all1, 2))
% 	n(itype) = size(oct_avg_all1, 2);
% 
% 	not_nan_ind = ~isnan(oct_avg_all);
% 	oct_avg = oct_avg_all(not_nan_ind);
% 
% 	if ndims(oct_mat) == 3
% 		oct_std_all = squeeze(std(oct_mat,0,[1 2], 'omitnan'))./sqrt(num_samples');
% 		oct_std = oct_std_all(not_nan_ind)';
% 	elseif ndims(oct_mat) == 2
% 		oct_std_all = std(oct_mat,0,1, 'omitnan')./sqrt(num_samples');
% 		oct_std = oct_std_all(not_nan_ind);
% 	end
% 	f_notnan = f(not_nan_ind);
% 
% 	hold on
% 	plot(f_notnan, oct_avg, 'color', 'k', 'LineWidth',1.5);
% 	patch([f_notnan flip(f_notnan)], [oct_avg-oct_std flip(oct_avg+oct_std)], ...
% 		colors{itype}, 'FaceAlpha',0.25, 'EdgeColor','none')
% 
% end
% yline(0)
% xline(0)
% xlim([-1.7 1.7])
% ylim([-1.5 2.5])
% for ii = 1:length(names)
% 	leg{2*ii-1} = '';
% 	leg{2*ii} = sprintf('%s, n=%d', names{ii}, n(ii));
% end
% 
% 
% CF_tests = -1:0.25:1;
% for ii = 1:9
% 	xline(CF_tests(ii), '--', 'color', [0.5 0.5 0.5])
% 	leg{length(names)*2+ii} = '';
% end
% set(gca,'XGrid','off','YGrid','on')
% hleg = legend(leg, 'Location','northeastoutside', 'FontSize',6);
% hleg.ItemTokenSize = [10,8];
% 
% end