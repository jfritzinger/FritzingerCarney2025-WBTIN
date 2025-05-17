function params = generate_TIN(params)
% TIN studies at desired center frequency. The bandwidth of the narrowband
% noise is the same with the critical bandwidth of the auditory filter.


%% Construct filters...
fn = params.Fs/2;
if params.bw_choice == 1
    lbw = params.CF * 2^(-1/6);
    hbw = params.CF * 2^(1/6);
elseif params.bw_choice ==2
    lbw = 100;
    hbw = 10e3;
elseif params.bw_choice == 3
    lbw = params.CF - 50;
    hbw = params.CF + 50;
elseif params.bw_choice == 4
    lbw = params.CF * 2^(-3/2); % JBF 6/13/22 spec_slide WB-TIN is 3 octave noise
    hbw = params.CF * 2^(3/2);
end
params.bw = hbw - lbw;
b_bp = fir1(5000,[lbw,hbw]/fn);         % Bandpass filter

%% Estimate time required for this DSID

nSNRs = length(params.SNR);
SNRs = params.SNR';
nNos = length(params.No);
nstim = (1 + nSNRs) * nNos;
ind1=ones(nSNRs+1,1)*params.No;
ind2=[-inf;SNRs]*ones(1,nNos);
spl_matrix = [ind1(:) ind2(:)];
params.spl_matrix = spl_matrix;
clear ind1 ind2;


%% Create the Stimulus Gating function
npts = floor(params.dur*params.Fs);
gate = tukeywin(npts,2*params.ramp_dur/params.dur); % Raised cosine ramps
t = (0:1:npts-1)/params.Fs;


%% Generate signal and noise
signal = sin(2*pi*params.CF*t');


%% Generate stimuli for all presentations
stimuli = zeros(npts,nstim*params.mnrep);
ipsi_stimuli = zeros(npts,nstim*params.mnrep);
presentation = 0;

rng('shuffle');
seed = randi(1e6,1); % to generate the "basic" seed for version 4. Was saving the seed structure.
for irep = 1:params.mnrep
    if params.choice~=3 && (irep ==1 || params.choice==2)
        params.seed(irep) = seed+irep;
        rng(params.seed(irep));
        noise = randn(npts,1);
        params.noise_first8_last8{irep} = [noise(1:8),noise(end-7:end)];
        noise = conv(noise,b_bp,'same');
    end

    stim_list = randperm(nstim);
    for istim = 1:nstim
        presentation = presentation + 1;

        if params.choice == 3
            params.seed(presentation) = seed+presentation;
            rng(params.seed(presentation));
            noise = randn(npts,1);
            params.noise_first8_last8{presentation} = [noise(1:8),noise(end-7:end)];
            noise = conv(noise,b_bp,'same');
        end

        irand = stim_list(istim);

        % Create matrices for No and SNR, the changing variables
        No = spl_matrix(irand,1);
        noise_spl = No + 10*log10(params.bw);
        SNR = spl_matrix(irand,2);

        % Scale noise in Pa
        noise_rms_Pa = (10^(noise_spl/20))*20e-6;   % scaler for desired sound level in Pa
        noise_Pa = noise.*(noise_rms_Pa/rms(noise));

        % Scales tone to correct SNR and level in Pa
        if SNR > -99
            signal_spl = No + SNR; % SNR is w.r.t. No
            signal_rms_Pa = (10^(signal_spl/20))*20e-6; % Note: signal is the tone
            signal_Pa = sqrt(2)* signal.*signal_rms_Pa; % first normalize then scale to desired tone level
        else
            signal_Pa = zeros(npts,1);
        end

        % Create contra and/or ipsi stimulus
        stim = noise_Pa + signal_Pa;
        params.stim(presentation,:) = stim.*gate;  % Gated stimulus with raised cosine ramps

        % Save
        params.mlist(presentation).SNR = SNR;
        params.mlist(presentation).No = No;
        params.mlist(presentation).noise_spl = noise_spl;
        params.mlist(presentation).irand = irand;
        params.mlist(presentation).rep = irep;
    end
end

end                       
