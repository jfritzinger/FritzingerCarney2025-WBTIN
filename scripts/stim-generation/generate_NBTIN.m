function [param] = generate_NBTIN(param)
% Generates the sliding wideband TIN stimulus for population model or
% single cell model
% J. Fritzinger, updated 6/13/2022
%
% Inputs: param.Fs: sampling rate
%         param.fpeak_mid: center frequency (tone frequency)
%         stp_oct: steps per octaves, 0 = population response to one
%         stimulus
%         bandwidth: noise bandwidth
%         spl: overall sound level
%         SNR: signal to noise ratio
%         Binmodes: 1 for contra, 2 for binaural


if param.physio == 1 % Copy from saved physio parameters 
	param.Fs = 48000;
	param.mnrep = param.nrep;

	if param.dur > 1
		param.dur = param.dur/1000; % stimulus duration in seconds.
	end

	rng(param.seed)
else
	lbw = param.fpeak_mid * 2^(-param.bandwidth/2);
	hbw = param.fpeak_mid * 2^(param.bandwidth/2);
	param.spl = param.No + 10*log10(hbw-lbw); % Overall level
	param.SNR = param.No + param.SNRNo - param.spl; % TIN SNR w.r.t. overall level

	% Generate noise for NB_noise
	rng('shuffle') % create a random seed
	rand_state = rng; % store it in the "seed field", a field of the random # generator
	param.seed = rand_state.Seed; % save this "basic" seed that will be used for waveform generation
	
end
NB_noise = flat_spec_noise(param.dur,param.Fs ); % NOTE: This noise will be frozen across reps and fpeaks!
disp(['Overall SPL Level = ' num2str(param.spl) 'dB']);

freq_lo = param.fpeak_mid/2; % range of  center-noise freqs: one octave below CF
freq_hi = param.fpeak_mid*2; %    "        "       "       : one octave above
if freq_hi > 20000 % limited to 20 kHz on upper end
	freq_hi = 20000;
end

% Calculate array of swept (center or edge) frequencies
if param.stp_otc == 0
	param.fpeaks = param.fpeak_mid;
else
	nfreqs = floor(param.stp_otc * log2(freq_hi/freq_lo))+1; % # of freqs, +1 makes array centered on CF
	param.fpeaks = freq_lo * 2.^((0:nfreqs-1)/param.stp_otc);
end
nstim = length(param.fpeaks)*length(param.SNR);

% Create the Stimulus Gating function
npts = floor(param.dur*param.Fs );
gate = tukeywin(npts,2*param.ramp_dur/param.dur); %raised cosine ramps

% Generate stimuli for all presentations
param.stim = zeros(param.mnrep*nstim, npts);
presentation = 0;
for irep = 1:param.mnrep
	%stim_list = randperm(nstim); % randomize stimuli for each rep
	stim_list = 1:nstim;

	% Create stimuli for each rep (irep) and each stimulus (istim)
	for istim = 1:nstim
		presentation = presentation + 1;

		% Set the RNG seed to stored "basic" value, and add a # to change it from rep to rep
		rng(param.seed + presentation); % is it true that we don't use random numbers after this point?

		% Compute one stimulus waveform.
		% npts = round(params.dur * fs);
		% magnitudes = ones(npts,1); % FLAT amplitude spectrum
		% phases = rand(npts,1);
		% noise = ifft(magnitudes + i*phases); % single noise is generated here - different bands will be presented for each fpeak
		% noise = real([0; noise(2:end)]); % !!! need to fix this! why is 1st point huge??
		% NOTE:  on 1/31/2020, moved this up above loop, so it will be frozen LHC
		% noise = flat_spec_noise(params.dur,fs); % NOTE: This noise will be different (different phases) for each rep & stimulus!

		nSNR = length(param.SNR);
		iSNR = mod(stim_list(istim)-1, nSNR)+1;
		SNR = param.SNR(iSNR);
		ifpeak = ceil((stim_list(istim))/nSNR);
		fpeak = param.fpeaks(ifpeak);
		contra_stim = NB_noise_gen(fpeak,NB_noise,param.bandwidth,...
			param.spl,param.Fs);
		if SNR ~= -inf  % add a tone at fpeak,  params.SNR ~= -inf was a bug prior to v.6
			tone = (20e-6 * 10.^((param.spl + SNR)/20) * sqrt(2) * ...
				sin(2*pi*fpeak*[1:length(NB_noise)]/param.Fs))';
			contra_stim = contra_stim + tone;
		end

		% Apply ramps 
		contra_stim = contra_stim.*gate;

		% Create stimulus matrix
		param.stim(presentation, :) = contra_stim;

		% SNR varies in WB/NB_noise, save SNR and iSNR
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

function noise = flat_spec_noise(dur,fs)
npts = round(dur*fs);
N = ceil(npts/2) - 1; % works for npts even or odd
phase = 2*pi*rand(N,1);
if rem(npts,2) % npts odd, No value at Nyquist freq.
	spec = [0;exp(1i*phase);flip(exp(-1i*phase))];
else % npts even, at Nyquist freq, mag = 1, and phase = -phase => phase = 0.
	spec = [0;exp(1i*phase);1;flip(exp(-1i*phase))];
end
noise = real(ifft(spec)); % real() is just insurance against tiny imaginary values.
end


function  contra_stim = NB_noise_gen(this_fpeak,noise,bandwidth, spl, fs)
% slide center frequency of a NB noise spectrum past CF, but phases will
% remain "fixed" on frequency axis. Create a single noise waveform above
%(with a tone at the CF), and just slide the NB filter past the noise to
%create each stimulus
fn = fs/2;
lbw = this_fpeak * 2^(-bandwidth/2);
hbw = this_fpeak * 2^(bandwidth/2);
b_bp = fir1(5000,[lbw,hbw]/fn);         % Bandpass filter
contra_stim = conv(noise,b_bp,'same');
scalar = (10^(spl/20))*20e-6;               % Target rms in Pa
contra_stim = contra_stim.*(scalar/rms(contra_stim));    % Noise at target level in Pa

end

