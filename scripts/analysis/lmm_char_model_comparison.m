%% lmm_char_model_comparison.m
% This script was used to compare linear mixed effects models, where the
% effects included characterizations of the neuron (MTF, CF group, etc). 
%
% J. Fritzinger
clear

%% Results, step-up (for each tone frequency

[base, datapath, savepath, ppi] = get_paths();
spreadsheet_name = 'Population_Data.xlsx';

analysisTable = readtable(fullfile(datapath, spreadsheet_name),...
	'PreserveVariableNames',true);
analysisTable.No = categorical(analysisTable.No);
analysisTable.No = reordercats(analysisTable.No, ["3", "23", "43"]);

%% Get rid of RM types with really small sample sizes & only 40 dB SNR

% Problems are with inhibitory, O, off, on, unusual, on/off
types = unique(analysisTable.RM);
for ii = 1:length(types)
	num_type(ii) = sum(strcmp(types{ii}, analysisTable.RM));
end
for ii = [2, 3, 4, 5, 6, 7]
	ind_inhibitory = strcmp(types{ii}, analysisTable.RM);
	analysisTable(ind_inhibitory,:) = [];
end

types = unique(analysisTable.SNR);
for ii = 1:2
	ind_inhibitory = types(ii)==analysisTable.SNR;
	analysisTable(ind_inhibitory,:) = [];
end

%%

% Only random effect
equation{1} = '~1+(1|Putative)';

% One fixed effect 
equation{2} = '~Tone_Freq+(1|Putative)'; % Best of this group
equation{3} = '~No+(1|Putative)';
equation{4} = '~Binmode+(1|Putative)';
equation{5} = '~RM+(1|Putative)';
equation{6} = '~MTF+(1|Putative)';
equation{7} = '~CF_Group+(1|Putative)';

% Three fixed effects
equation{8} = '~Tone_Freq*(No*RM)+(1|Putative)'; % Best of this group
equation{9} = '~Tone_Freq*(No*MTF)+(1|Putative)';
equation{10} = '~Tone_Freq*(No*CF_Group)+(1|Putative)'; 

equation{27} = '~Tone_Freq*(No*RM*CF_Group)+(1|Putative)';
equation{28} = '~Tone_Freq*(No*RM*CF_Group)+Tone_Freq*(No*MTF)+(1|Putative)';
equation{29} = '~Tone_Freq*(No*RM*CF_Group)+Tone_Freq*(No*MTF*CF_Group)+(1|Putative)';


equation{11} = '~Tone_Freq*(No*RM*MTF)+(1|Putative)';
equation{12} = '~Tone_Freq*(No*MTF*CF_Group)+RM+(1|Putative)'; % Best
equation{13} = '~Tone_Freq*(No*CF_Group*RM)+(1|Putative)';
equation{14} = '~Tone_Freq*(No*CF_Group*RM*MTF)+(1|Putative)';
equation{15} = '~Tone_Freq*(CF_Group*RM*MTF)+(1|Putative)';


equation{16} = '~Tone_Freq*(No*MTF*CF_Group)+(1|Putative)';
equation{17} = '~Tone_Freq*(No*MTF*CF_Group)+RM+(1|Putative)';
equation{18} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM+(1|Putative)'; % Best
equation{19} = '~Tone_Freq*No+Tone_Freq*MTF+Tone_Freq*CF_Group+Tone_Freq*RM+(1|Putative)';
equation{20} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM+Binmode+(1|Putative)';
equation{21} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM*CF_Group+(1|Putative)';
equation{22} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM*MTF+(1|Putative)';
equation{23} = '~Tone_Freq*No*MTF*CF_Group*RM+(1|Putative)'; % Error
equation{24} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM*MTF+Tone_Freq*No*RM*CF_Group+(1|Putative)'; % Best if cutting down to 2 RM groups
equation{25} = '~Tone_Freq*No*MTF*CF_Group+Tone_Freq*No*RM*MTF+Tone_Freq*No*RM*CF_Group+Tone_Freq*RM*CF_Group*MTF+(1|Putative)'; % Best if cutting down to 2 RM groups
equation{26} = '~Tone_Freq*(RM*MTF*CF_Group)+(1|Putative)'; % Best if cutting down to 2 RM groups


%%

a = 'Tone_Freq*No + Tone_Freq*CF_Group + Tone_Freq*SNR + Tone_Freq*MTF + Tone_Freq*RM';
b = '+Tone_Freq*No*CF_Group + Tone_Freq*SNR*CF_Group';
c = '+Tone_Freq*No*RM + Tone_Freq*SNR*RM';
d = '+Tone_Freq*No*MTF + Tone_Freq*SNR*MTF';
e = '+Tone_Freq*No*SNR*MTF + Tone_Freq*No*SNR*CF_Group + Tone_Freq*No*SNR*RM';
f = '+Tone_Freq*RM*MTF';

randeffect = '+ (1 | Putative)';

%% Test 
ref =  11; % 9, 13, 16
test = 24; 

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

full_equation = ['Rate' equation{24}];
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