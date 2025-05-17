function [dog, gauss_exc, gauss_inh] = createDoG(DOGparams, fpeaks)
%CREATEDOG Generate Difference of Gaussians (DoG) model for neural responses
%   [dog, gauss_exc, gauss_inh] = createDoG(DOGparams, fpeaks) creates a
%   Difference of Gaussians
%
%   Inputs:
%       DOGparams - 6-element vector of model parameters:
%           [s_exc, s_inh, sigma_exc, sigma_inh, CF_exc, CF_inh]
%           s_exc:    Excitation strength (a.u.)
%           s_inh:    Inhibition strength (a.u.)
%           sigma_exc: Excitation bandwidth (Hz)
%           sigma_inh: Inhibition bandwidth (Hz)
%           CF_exc:    Excitation center frequency (Hz)
%           CF_inh:    Inhibition center frequency (Hz)
%       fpeaks - Frequency values to evaluate model (vector)
%
%   Outputs:
%       dog        - Difference of Gaussians model (s_exc*gauss_exc - s_inh*gauss_inh)
%       gauss_exc  - Normalized excitation Gaussian component
%       gauss_inh  - Normalized inhibition Gaussian component


% Unpack parameters
s_exc = DOGparams(1);
s_inh = DOGparams(2);
sigma_exc = DOGparams(3);
sigma_inh = DOGparams(4);
CF_exc = DOGparams(5);
CF_inh = DOGparams(6);

% Create base Gaussian components
gauss_exc = normpdf(fpeaks, CF_exc, sigma_exc);
gauss_inh = normpdf(fpeaks, CF_inh, sigma_inh);

% Normalize and scale components
gauss_exc = s_exc*(gauss_exc./max(gauss_exc));
gauss_inh = s_inh*(gauss_inh./max(gauss_inh));

% Calculate DoG response
dog = gauss_exc - gauss_inh;

end

