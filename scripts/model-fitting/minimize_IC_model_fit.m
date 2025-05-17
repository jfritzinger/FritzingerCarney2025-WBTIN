function mse = minimize_IC_model_fit(data, AN, stim_params, model_params, x)
% minimize_IT_model_fit.m runs the lateral SFIE variable BMF function and
% then calculated the MSE between the model output and data for MTF, TIN,
% and WB-TIN. 
%
% J. Fritzinger

% Run model 
CS_params = [x(1:2) 0.001];
BMFs = x(3:5);
nstim = size(stim_params, 2);
model_outputs = cell(nstim, 1);

parfor istim = 1:nstim
	param = stim_params{istim};
	an_sout = squeeze(AN{istim}.an_sout);
	an_sout_lo = squeeze(AN{istim}.an_sout_lo);
	an_sout_hi = squeeze(AN{istim}.an_sout_hi);
	model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
end

% Analyze model output
for istim = 1:nstim
	param = stim_params{istim};
	model_output = model_outputs{istim};
	switch param.type
		case 'typMTFN'
			[~, model_MTF, ~, ~] = plotMTF(param, model_output.avIC, 0);
		case 'TIN'
			[~, model_TIN,~] = plotTIN(param, model_output.avIC, 0);
		case 'SPEC_slide'
			[~, model_WB, ~] = plotWBTIN(param, model_output.avIC, 0);
	end
end

%% Calculating MSE 

% Normalize to match real data 
data_MTF = data(1:26);
data_TIN = data(27:29);
data_WB = data(30:end);

% Minimize -1*R
r_MTF = corrcoef(model_MTF, data_MTF);
r_TIN = corrcoef(model_TIN, data_TIN);
r_WB = corrcoef(model_WB, data_WB);
r_mean = (r_MTF(1,2) + r_TIN(1,2) + r_WB(1,2))/3;
mse = -1* r_mean;

%% Display current MSE and parameters in command window

r_MTF = corrcoef(model_MTF, data_MTF);
r_TIN = corrcoef(model_TIN, data_TIN);
r_WB = corrcoef(model_WB, data_WB);
fprintf('S1=%0.03f, S2=%0.03f, D=1ms, BMFs=[%0.03f %0.03f %0.03f]\n',...
	x(1), x(2), x(3), x(4), x(5)) 
fprintf('R=%0.2f, RMTF=%0.2f, RTIN=%0.2f, RWB=%0.2f\n',...
	mse, r_MTF(1,2), r_TIN(1,2), r_WB(1,2))

end
