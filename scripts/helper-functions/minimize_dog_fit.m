function mse = minimize_dog_fit(rate, fpeaks, DOGparams, mse_type)
% Calculate error metric for Difference of Gaussians (DoG) fit
%   mse = MINIMIZE_DOG_FIT(rate, fpeaks, DOGparams, mse_type) computes the
%   error between observed neural response rates and a DoG model prediction.
%
%   Inputs:
%       rate      - Observed response rates (vector)
%       fpeaks    - Frequency values corresponding to rate measurements (vector)
%       DOGparams - Parameters for Difference of Gaussians model (struct/vector)
%       mse_type  - Error calculation method (1 or 2)
%                  1: Normalized L2-norm (Euclidean distance)
%                  2: Standard Mean Squared Error (MSE)
%
%   Output:
%       mse - Calculated error metric (scalar)

% Calculate DoG 
[dog, ~, ~] = createDoG(DOGparams, fpeaks);

if mse_type == 1 % Minimize squared error of Euclidian distance (L2)
	mse = norm(rate-dog,2)/norm(rate);

elseif mse_type == 2 % Compute the mean of the squared errors
	mse = 1/length(rate) * sum((rate - dog).^2);

end

end