function data  = analyzeRM(params)
% Analyze response maps, pure tone responses at different frequencies and
% sound levels. 
%   data = analyzeRM(params) processes neural responses to characterize
%   response map and outputs average rates to plot.
%
%   Inputs:
%       params - Structure containing experimental parameters for RM
%
%   Outputs:
%       data - Structure containing response map analysis results
%           Fields:
%           .rates  - Mean firing rate matrix (frequency × SPL) in spikes/sec
%           .spont  - Spontaneous rate (mean response at lowest SPL)
%           .freqs  - Unique frequencies analyzed (Hz)
%           .SPL    - Unique sound pressure levels analyzed (dB SPL)
%
%   See also: accumstats

this_ds = params.stims.dsid==params.dsid;

[spls,~,si] = unique(double([params.list.spl]).');
num_spls = length(spls);
[freqs,~,fi] = unique(double([params.list.freq]).');
num_freqs = length(freqs);

rate_size = [num_freqs,num_spls];
num_spikes = params.cluster.num_spikes_peri(this_ds);
spike_rates = num_spikes*1000/params.dur; % spikes/sec

if length(si) == length(spike_rates)
    [rates,~,~,~] = accumstats({fi,si},spike_rates, rate_size);
	spont = mean(rates(:,1));  
end

data.rates = rates;
data.spont = spont;
data.freqs = freqs;
data.SPL = spls;

end