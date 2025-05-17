function [param] = generate_WBTIN(param)
% Generates the sliding wideband TIN stimulus for population model or
% single cell model
% J. Fritzinger, updated 6/13/2022
%
% Inputs: fs: sampling rate
%         Fc: center frequency (tone frequency)
%         stp_oct: steps per octaves, 0 = population response to one
%         stimulus
%         bandwidth: noise bandwidth
%         spl: overall sound level
%         SNR: signal to noise ratio
%         Binmodes: 1 for contra, 2 for binaural


% Parameters
if param.physio == 1 % Copy from saved physio parameters 
	param.mnrep = param.nrep;
    param.CF = param.fpeak_mid;
    param.stp_oct = param.stp_otc;
	param.Fs = 48000;
	WB_seed = param.WB_seed;
	if param.dur > 1
		param.dur = param.dur/1000; % stimulus duration in seconds.
	end
	seed = param.seed;
else
	% Noise Seed
    rng('shuffle') % create a random seed
    WB_rand_state = rng; % store it in the "seed field", a field of the random # generator
    WB_seed = WB_rand_state.Seed; % save this "basic" seed that will be used for waveform generation

	% Random presentation seed
	rng('shuffle') % create a random seed
	rand_state = rng; % store it in the "seed field", a field of the random # generator
	seed = rand_state.Seed; % save this "basic" seed that will be used for waveform generation
	if param.dur > 1
		param.dur = param.dur/1000; % stimulus duration in seconds.
	end
end

if param.version >= 5
    param.freq_lo = param.fpeak_mid*2^(-param.range/2); % range of  center-noise freqs
    param.freq_hi = param.fpeak_mid*2^(param.range/2); %    "        "       "     
else
    param.freq_lo = param.CF/2; % range of  center-noise freqs: one octave below CF
    param.freq_hi = param.CF*2; %    "        "       "       : one octave above
end
if param.freq_hi > 20000 % limited to 20 kHz on upper end
    param.freq_hi = 20000;
end

% Calculate stimulus frequencies
if param.stp_otc == 0
    fpeaks = [0 param.CF];
	param.fpeaks = fpeaks;
elseif param.physio == 1 % Copy from saved physio parameters 
	fpeaks = param.fpeaks;
else
    nfreqs = floor(param.stp_otc * log2(param.freq_hi/param.freq_lo))+1; % # of center-noise freqs
    fpeaks = [0 param.freq_lo * 2.^((0:nfreqs-1)/param.stp_otc)];
	param.fpeaks = fpeaks;
end
nfreqs = length(fpeaks);

% Calculate stimulus SNRs
nSNRs = length(param.SNR);

% Estimate time required
nstim = nfreqs*nSNRs; % just for display
npts = floor(param.dur*param.Fs);
param.steps = nstim;

%% Generate stimuli for all presentations

% Display
lbw = param.CF * 2^(-param.bandwidth/2);
hbw = param.CF * 2^(param.bandwidth/2);
overall_level = param.No + 10*log10(hbw-lbw);
%disp(['Overall SPL Level = ' num2str(overall_level) 'dB']);
%assert(overall_level<85, 'Overall noise level is above 85dB!');

presentation = 0; %this value is used as an index for storing a stumulus presentation in the 3rd dimenstion of 'stimuli'
param.stim = zeros(param.mnrep*nstim, npts);
rng(seed) % create a random seed
for irep = 1:param.mnrep
	if param.stp_otc == 0 % not randomized if population model
		stim_list = 1:nstim;
	elseif param.physio == 1 % Copy from saved physio parameters
		stim_list = [param.list.irand];
	else
		stim_list = randperm(nstim);
	end
    WB_noise = WB_noise_gen(param.fpeak_mid, param.dur, ...
                 WB_seed+irep, param.bandwidth, overall_level, param.Fs);  
    param.overall_level = overall_level;
    for istim = 1:nstim
        presentation = presentation + 1;
        
        % Set the RNG seed to stored "basic" value, and add a # to change it from rep to rep
        rng(seed + presentation);
        
        % Get the peak (or center, or edge) freq for this presentation.
        iSNR = mod(stim_list(istim)-1,nSNRs)+1;   % the "-1" makes sure that list of SPLS starts with low value, and "+1" starts values at 1 instead of 0
        SNR = param.SNR(iSNR);
        ifpeak = ceil((stim_list(istim))/nSNRs);
        fpeak = fpeaks(ifpeak);
        
        % Add tone
        if fpeak ~= 0
            tone = (20e-6 * 10.^((param.No + SNR)/20) * sqrt(2) * ...
                sin(2*pi*fpeak*[1:length(WB_noise)]/param.Fs))';
            stim = WB_noise + tone;
        else
            stim = WB_noise;
        end
        
        % Apply ramps.
        gate = tukeywin(npts,2*param.ramp_dur/param.dur);
        param.stim(presentation, :) = stim.*gate;
        
        % Save freq and SNR lists
        param.mlist(presentation).fpeak = fpeak;
        param.mlist(presentation).ifpeak = ifpeak;
        param.mlist(presentation).SNR = SNR;
        param.mlist(presentation).iSNR = iSNR;
        param.mlist(presentation).rep = irep;
        param.mlist(presentation).irand = stim_list(istim);
        
    end
end

end

%% %%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  contra_stim = WB_noise_gen(CF,dur, WB_seed, bandwidth, No_overall, fs)
% Creates a wideband flat spectrum noise and filters the noise to the proper bandwidth
% J. Fritzinger, updated 2/17/2022

% Creates noise 
npts = round(dur*fs);
N = ceil(npts/2) - 1; % works for npts even or odd
rng(WB_seed);
phase = 2*pi*rand(N,1);
if rem(npts,2) % npts odd    
    spec = [0;exp(1i*phase);flip(exp(-1i*phase))]; % No value at Nyquist freq.
else % npts even   
    spec = [0;exp(1i*phase);1;flip(exp(-1i*phase))]; % At Nyquist freq, mag = 1, and phase = -phase => phase = 0.
end
noise = real(ifft(spec)); % real() is just insurance against tiny imaginary values.

% Filters noise 
fn = fs/2;
lbw = CF * 2^(-bandwidth/2);
hbw = CF * 2^(bandwidth/2);

if hbw > 20000 % limited to 20 kHz on upper end
    hbw = 20000;
end

b_bp = fir1(5000,[lbw,hbw]/fn);         % Bandpass filter
contra_stim = conv(noise,b_bp,'same');
scalar = (10^(No_overall/20))*20e-6;               % Target rms in Pa
contra_stim = contra_stim.*(scalar/rms(contra_stim));    % Noise at target level in Pa

end
