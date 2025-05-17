function e = calculateEllipseSTD(X)
%Calculate standard deviation ellipse for 2D data
%   e = calculateEllipseSTD(X) computes the standard deviation ellipse 
%   (covariance ellipse) for a set of 2D data points.
%
%   Inputs:
%       X - N×2 matrix of 2D data points (N samples)
%
%   Outputs:
%       e - 2×100 matrix containing ellipse coordinates (x; y)

% Subtract mean
num_samples = length(X);
Mu = mean(X);
X0 = bsxfun(@minus, X, Mu);

% Covariance matrix via eigenvalue decomposition
[V, D] = eig(X0'*X0 ./ (num_samples-1));  % Covariance matrix calculation
[D, order] = sort(diag(D), 'descend');    % Sort eigenvalues/vectors
D = diag(D);                               % Maintain matrix form
V = V(:, order);                           % Reorder eigenvectors

% Generate unit circle and transform
t = linspace(0, 2*pi, 100);               % 100-point circle parameterization
e = [cos(t); sin(t)];                     % Unit circle coordinates
VV = V*sqrt(D);                           % Scale eigenvectors by std devs
e = bsxfun(@plus, VV*e, Mu');             % Rotate and translate to data space

end
