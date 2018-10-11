% Generate Sobol LP sequence for sensitivity analysis of ebullition model
% 1 Ks:        thermal conductivity of soil constituent (0.25~2.9 W/(m*K))
% 2 Cps:       heat capacity of soil constituent (750~1930 J/(kg*K))
% 3 Porosity:  porosity of lake sediment (0.3~0.6)
% 4 Rous:      density of soil solid particle (1500~2700 kg/m3)
% 5 OQ10:      CH4 oxidation Q10 (1.4~3.5)
% 6 Qch4:      Oxidation potential when substrates are not limited (0.1~100 umol/(m3*s))
% 7 Kch4:      Michaelis-Menten constant (1e3~6.62e4 umol/m3)
% 8 Ko2:       Michaelis-Menten constant (1e3~2.0e5 umol/m3)
% 9 Re:        ebullition rate (1.85e-4~9.25e-4 s-1)
%10 PQ10n:     Methane production Q10 of acetate fermentation (1.7~16)
%11 Rcn:       the fraction of decomposed recalcitrant carbon (6.342e-11~6.342e-10 s-1)
%12 DMP:       recalcitrant organic matter dampening rate (1.0~10.0 m-1)
%13 Rco:       the fraction of decomposed labile carbon (6.342e-11~3.171e-9 s-1)
%14 PQ10o:     Methane production Q10 of CO2 reduction (1~3.5)
%15 Rca:       the fraction of aerobic decomposed carbon (2.3e-10~5.8e-9 s-1)
%16 Vchs:      Chla-specific light saturated growth rate (18.26~200 mg C mg Chl-1 d-1)
%17 Vchl:      Chla-specific light saturated growth rate (18.26~200 mg C mg Chl-1 d-1)
%18 Klrs:       metabolic loss rate coefficient (0.04~0.125 day-1)
%19 Klrl:       metabolic loss rate coefficient (0.04~0.125 day-1)
%20 RDOMaq:    aquatic DOM microbial degradation rate (0.01~0.1 d-1)
%21 RDOMtr:    terrestrail DOM microbial degradation rate (5e-4~0.05 d-1)
%22 RDOMno:    aquatic DOM anearobic degradation rate (1e-5~2.5e-3 d-1)
%23 DOCwt:     active layer DOC concentration (0.66~7.78 mol/m3)
%24 Ero:       thermokarst erosion multiplier (0.02~1.81 m/yr)
%25 Wro:       wind erosion multiplier (0.005~0.015 m/yr/W)
%26 phAlphas:  initial slope of P-E curve (1.5e3~5.6e4 (mol photons m-2 s-1)-1 d-1)
%27 phAlphal:  initial slope of P-E curve (2e3~3.6e4 (mol photons m-2 s-1)-1 d-1)
%28 phBetas:   photoinhibition parameter (2e2~2e3 (mol photons m-2 s-1)-1 d-1)
%29 phBetal:   photoinhibition parameter (2e2~2e3 (mol photons m-2 s-1)-1 d-1)
%30 Ksrps:     half-saturation for phosphorus limitation (3.226e1~4.839e3 umol/m3)
%31 Ksrpl:     half-saturation for phosphorus limitation (3.226e1~3.226e2 umol/m3)
%32 Feta:      light attenuation correction factor for chla (0.1~10)

optpnt = [2.830968e+00, 1.769980e+03, 3.247692e-01, 2.240472e+03, ...
   1.429714e+00, 4.806569e-01, 1.408063e+04, 7.810067e+04, 8.951114e-04, ...
   3.941954e+00, 1.490403e-10, 3.528460e+00, 2.149145e-10, 1.006950e+00, ...
   5.791722e-09, 1.954113e+02, 5.392818e+01, 4.844265e-02, 4.844265e-02, ...
   9.610503e-02, 2.800566e-03, 7.276936e-03, 1.026985e+01, 8.285088e-02, ...
   1.000000e-02, 5.178555e+04, 2.970566e+04, 3.885597e+02, 4.462115e+02, ...
   1.580650e+02, 3.612900e+02, 3.708900e+00];

test_indice = [3, 11, 12, 16, 18, 23, 32];

%test_indice = [3, 16, 18, 20, 23];

%parmin = [6.342e-11, 1.000e+00, 1.826e+01, 1.826e+01, 1.000e-01];
 
%parmax = [6.342e-10, 1.000e+01, 2.000e+02, 2.000e+02, 1.000e+01];
%parmin = [3.000e-01, 1.836e+01, 4.000e-02, 1.000e-02, 6.600e-01];
%parmax = [6.000e-01, 2.000e+02, 1.250e-01, 1.000e-01, 2.000e+01];

parmin = [3.000e-01, 6.342e-11, 1.000e+00, 1.836e+01, 4.000e-02, 6.600e-01, 1.000e-01];
parmax = [6.000e-01, 6.342e-10, 1.000e+01, 2.000e+02, 1.250e-01, 2.000e+01, 1.000e+01];

nparam = length(optpnt);
ntest = length(test_indice);

nsample = 10000;
nskip = 2^min(ntest,15);
sample_1 = sobol_dataset( ntest, nsample, nskip );
sample_2 = sobol_dataset( ntest, nsample, 2*nskip );
sample_1 = sample_1';
sample_2 = sample_2';
for ii = 1 : ntest
   sample_1(:,ii) = log(parmin(ii))*ones(nsample,1) + ...
      (log(parmax(ii)) - log(parmin(ii)))*sample_1(:,ii);
   sample_1(:,ii) = exp(sample_1(:,ii));
   sample_2(:,ii) = log(parmin(ii))*ones(nsample,1) + ...
      (log(parmax(ii)) - log(parmin(ii)))*sample_2(:,ii);
   sample_2(:,ii) = exp(sample_2(:,ii));
end

% Extended Sobol method needs to construct nparam + 2 samples including sample_1 and
% sample_2 (A. Saltelli (2002): Making best use of model evaluations to compute 
% sensitivity indices, Computer Physics Communications, 145, 280-297).
fmtstr1 = repmat(' %14.6e', 1, nparam);
fmtstr = ['%5.0d',fmtstr1,'\n'];
% save samples to model_parameters.dat
fid = fopen('./model_parameters.dat','w');
for ii = 1 : nsample
   sampleline = optpnt;
   sampleline(test_indice) = sample_2(ii,:);
   fprintf(fid, fmtstr, ii, sampleline);
end
fclose(fid);
