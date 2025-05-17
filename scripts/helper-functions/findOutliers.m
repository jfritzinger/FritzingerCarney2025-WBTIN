function outliers = findOutliers(data)

data2 = data(:,1:2);

% Calculate the Mahalanobis distance for each data point
mahal_dist = mahal(data2(:,1), data2(:,2));

% Set a threshold for outlier detection (e.g., 95% confidence interval)
alpha = 0.05;
chi_sq_critical = chi2inv(1 - alpha, 2); % For 2 dimensions

% Identify outliers based on the threshold
outliers = find(mahal_dist > chi_sq_critical);

end