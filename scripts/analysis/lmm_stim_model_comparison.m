%% lmm_char_model_comparison.m
% This script was used to compare linear mixed effects models where the
% effects were limited to stimulus parameters (SNR, No, presentation ear,
% etc). 
%
% J. Fritzinger
clear

%% All: Run 0 vs 1 effect  

[base, datapath, savepath, ppi] = get_paths();
spreadsheet_name = 'Population_Data.xlsx';

analysisTable = readtable(fullfile(datapath, spreadsheet_name),...
	'PreserveVariableNames',true);
analysisTable.No = categorical(analysisTable.No);
analysisTable.No = reordercats(analysisTable.No, ["3", "23", "43"]);
analysisTable.SNR = categorical(analysisTable.SNR);

% Only random effect
equation{1} = '~1+(1|Putative)';
equation{2} = '~Tone_Freq+(1|Putative)';

% One fixed effect 
equation{3} = '~Tone_Freq*No+(1|Putative)';
equation{4} = '~Tone_Freq*SNR+(1|Putative)';
equation{5} = '~Tone_Freq*Binmode+(1|Putative)';

% Two fixed effects, added
equation{6} = '~Tone_Freq*(No+SNR)+(1|Putative)';
equation{7} = '~Tone_Freq*(No+Binmode)+(1|Putative)';
equation{8} = '~Tone_Freq*(Binmode+SNR)+(1|Putative)';

% Two fixed effects, multiplied
equation{9} = '~Tone_Freq*(No*SNR)+(1|Putative)';
equation{10} = '~Tone_Freq*(No*Binmode)+(1|Putative)';
equation{11} = '~Tone_Freq*(Binmode*SNR)+(1|Putative)';

% Three fixed effects
equation{12} = '~Tone_Freq*(No+SNR+Binmode)+(1|Putative)';
equation{13} = '~Tone_Freq*(No*SNR+Binmode)+(1|Putative)';
equation{14} = '~Tone_Freq*(No+SNR*Binmode)+(1|Putative)';
equation{15} = '~Tone_Freq*(No*Binmode+SNR)+(1|Putative)';
equation{16} = '~Tone_Freq*(No*SNR*Binmode)+(1|Putative)';

%% Test 
ref =  13; % 9, 13, 16
test = 16; 

full_equation = ['Rate' equation{ref}];
mdl = fitlme(analysisTable, full_equation,'FitMethod','ML', 'DummyVarCoding','effects');
R2 = mdl.Rsquared.Ordinary;

full_equation = ['Rate' equation{test}];
altmdl = fitlme(analysisTable, full_equation,'FitMethod','ML', 'DummyVarCoding','effects');
altR2 = altmdl.Rsquared.Ordinary;

disp('----------------------------------------------------------------------------------------------------------')
a = compare(mdl, altmdl);
disp(a)
fprintf('Model 1: R^2 = %0.6f\nModel 2: R^2 = %0.6f\n', R2, altR2)
disp('----------------------------------------------------------------------------------------------------------')

%% Example 

full_equation = ['Rate' equation{16}];
mdl = fitlme(analysisTable, full_equation,'FitMethod','REML', 'DummyVarCoding','effects');
R2 = mdl.Rsquared.Ordinary;
disp(mdl)
anova(mdl, 'DFMethod','satterthwaite')

figure
plotResiduals(mdl, 'fitted')
figure
plotResiduals(mdl,'probability')

F = fitted(mdl);
R = response(mdl);
figure();
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')

