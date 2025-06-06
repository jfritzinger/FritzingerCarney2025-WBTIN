function [ic_sout_BE,ic_sout_BS,ic_sout_hybrid, cn_sout] = SFIE_Hybrid_BMF(an_sout, BMF, fs)
% Expanded original SFIE model (Nelson & Carney, 2004 JASA) to include
% Band-Suppressed MTF (i.e. Low-Pass/Notch or High-pass) by adding a cell 
% that is excited by the CN input and inhibited by the Band-Enhanced Cell 
% (i.e. Bandpass cell) - see eNeuro, Carney et al., 2015.
% Adjustable BMF, per Ken's model
%
% Added parameters for hybrid cell, per 2016 model. J. Fritzinger. 
%
% Model from 2016 paper 

%% Model Parameters 

% CN MODEL PARAMETERS:
tau_ex_cn = 0.5e-3;         % CN exc time constant
tau_inh_cn = 2.0e-3;         % CN inh time constant
cn_delay = 1.0e-3;        % "disynaptic inhibition delay" (all ANFs excitatory)
inh_str_cn = 0.6;       % re: excitatory strength == 1
afamp_cn = 1.5;         % alpha function area --> changes RATE of output cell

% IC MODEL PARAMETERS:
% BMF-dependent SFIE parameters
tau_ex_ic = 1/(10*BMF); %[0.00025 0.0005 0.001 0.002];		% Time constant excitation in seconds
tau_inh_ic = tau_ex_ic*1.5;								% Time constant inhibition in seconds
ic_delay_inh = tau_ex_ic*2;% Delay of inhibition in seconds
afamp_ic = 1;             % alpha function area --> changes RATE of output IC BE cell 
inh_str_ic = 0.9; % inhibitory strength

% BS parameters
inh_str_bs = 4;  
tau_inh_bs = tau_inh_ic; %1.0e-3; % relatively long inhibition, from BE to BS
ic_delay_bs = 1.0e-3;  % Delay from BE to BS cell (local)
Aex = 0.5; %0.3; % Rate Scalar for BS cell; note that this is effectively multiplied by afamp_ic (for Table in eNeuro)
 
% Hybrid parameters 
tau_inh_hybrid = 6.0e-3;
AinhBE = 6.0;
AinhBS = 0.7;
ic_delay_hybrid = 1.5e-3;
Aex_hybrid = 1.5;

% tau_inh_hybrid = 6.0e-3;
% AinhBE = 3;
% AinhBS = 3;
% ic_delay_hybrid = 1.5e-3;
% Aex_hybrid = 1.5;

%% CN model:
% Generate frequency-domain equivalent of alpha functions
[B1, A1] = get_alpha_norm(tau_ex_cn, fs, 1);
[B2, A2] = get_alpha_norm(tau_inh_cn, fs, 1);
cn_ex = [afamp_cn*(1/fs)*(filter(B1, A1, an_sout)) zeros(1,fs*cn_delay)];
cn_inh = [zeros(1,fs*cn_delay) afamp_cn*inh_str_cn*(1/fs)*(filter(B2, A2,an_sout))];

% final CN model response:
cn_sout = ((cn_ex-cn_inh) + abs(cn_ex-cn_inh))/2;   % subtract inhibition from excitation and half-wave-rectify
%cn_t = 0:(length(cn_sout)-1)/fs;        % time vector for plotting CN responses

%% IC Model #1: (SFIE; Bandpass MRF)
% Generate alpha functions for BP IC model (same as CN model, but with different taus) See Nelson & Carney 2004
[B3, A3] = get_alpha_norm(tau_ex_ic, fs, 1);
[B4, A4] = get_alpha_norm(tau_inh_ic, fs, 1);
ic_lp_ex1 = [afamp_ic*(1/fs)*(filter(B3, A3, cn_sout)) zeros(1,floor(fs*ic_delay_inh)) zeros(1,floor(fs*ic_delay_bs))];
ic_lp_inh1 = [zeros(1,floor(fs*ic_delay_inh)) afamp_ic*inh_str_ic*(1/fs)*(filter(B4, A4, cn_sout)) zeros(1,floor(fs*ic_delay_bs))];
ic_sout_BE = ((ic_lp_ex1-ic_lp_inh1) + abs(ic_lp_ex1-ic_lp_inh1))/2; % half-wave rectified; standard SFIE model

%%  Band-suppressed cell (see Carney et al., 2015)
[B5, A5] = get_alpha_norm(tau_inh_bs, fs, 1);
ic_bs_ex = Aex * ic_lp_ex1; % add zeros at end to match lengths  
ic_bs_inh = Aex*inh_str_bs*(1/fs)*(filter(B5, A5,ic_sout_BE));
ic_sout_BS = ((ic_bs_ex-ic_bs_inh) + abs(ic_bs_ex-ic_bs_inh))/2; % half-wave rectified

%% Hybrid MTF 
[B7, A7] = get_alpha_norm(tau_inh_hybrid, fs, 1);
ic_hybrid_ex = Aex_hybrid*[ic_lp_ex1 zeros(1,floor(fs*ic_delay_hybrid))]; % add zeros at end to match lengths  
ic_hybrid_inhBE = Aex_hybrid*[zeros(1,floor(fs*ic_delay_hybrid)) AinhBE*(1/fs)*(filter(B7, A7, ic_sout_BE))];
ic_hybrid_inhBS = Aex_hybrid*[zeros(1,floor(fs*ic_delay_hybrid)) AinhBS*(1/fs)*(filter(B7, A7,ic_sout_BS))];
ic_sout_hybrid = ((ic_hybrid_ex-ic_hybrid_inhBS-ic_hybrid_inhBE) + abs(ic_hybrid_ex-ic_hybrid_inhBS-ic_hybrid_inhBE))/2; % half-wave rectified

% figure;
% plot(ic_lp_ex1)
% hold on
% plot(ic_lp_inh1)
% plot(ic_lp_ex1-ic_lp_inh1)
% yline(0, 'k')
% %plot(ic_sout_BE)
% title('BE Components')
% legend('Excitatory', 'Inhibitory', 'Combined')

% subplot(1, 2, 2);
% plot(ic_bs_ex)
% hold on
% plot(ic_bs_inh)
% plot(ic_bs_ex-ic_bs_inh)
% %plot(ic_sout_BS)
% title('BS Components')
% legend('Excitatory', 'Inhibitory', 'Combined')


end