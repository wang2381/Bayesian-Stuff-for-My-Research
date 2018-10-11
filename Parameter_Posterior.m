function [ optimum ] = Parameter_Posterior( paramId, obsIds )
% calculate the posterior probability of a parameter by bayesian inference
% and return the CDF of posterior distribution and the optimum parameter

%parmin = [2.500e-01, 7.500e+02, 3.000e-01, 1.500e+03, 1.400e+00, ...
%   1.000e-01, 1.000e+03, 1.000e+03, 1.850e-04, 1.700e+00, ...
%   6.342e-11, 1.000e+00, 6.342e-11, 1.000e+00, 1.000e+00, ...
%   1.000e+00];

%parmax = [2.900e+00, 1.930e+03, 6.000e-01, 2.700e+03, 3.500e+00, ...
%   1.000e+02, 6.620e+04, 2.000e+05, 9.250e-04, 1.600e+01, ...
%   6.342e-10, 1.000e+01, 3.171e-09, 3.500e+00, 8.000e+00, ...
%   1.000e+01];

parmin = [3.000e-01, 6.342e-11, 1.000e+00, 1.836e+01, 4.000e-02, 6.600e-01, 1.000e-01];
parmax = [6.000e-01, 6.342e-10, 1.000e+01, 2.000e+02, 1.250e-01, 2.000e+01, 1.000e+01];

%dir = '/scratch/rice/w/wang2381/paper4/';
dir = './';
nparam = length(parmin);

nsample = 2000;
hi = zeros(nsample,1);
filename = strcat(dir, 'SA', '.nc');
for ii = 1 : length(obsIds)
%   filename = strcat(dir, 'SA', num2str(obsIds(ii)), '.nc');
   ncid = netcdf.open(filename, 'NC_NOWRITE');
   varid = netcdf.inqVarID(ncid, 'sa');
   filval = netcdf.getAtt(ncid,varid,'_FillValue');
   sa_arr = netcdf.getVar(ncid, varid);
   sa_arr = normalize(sa_arr(ii, 1:end));
   sa_arr(sa_arr==filval) = -9999999;
   %
   sa_arr = reshape(sa_arr, [nsample,1]);
   hi = hi + sa_arr;
end
netcdf.close(ncid);
mhi = max(hi);
hi = exp(hi - mhi);
% construct resampling cumulative probability distribution
hisum = sum(hi) * ones(nsample,1);
rhi = hi ./ (hisum - hi);
mrhi = sum(rhi);
rcdf = zeros(nsample,1);
for ii = 1 : nsample
   rcdf(ii) = sum(rhi(1:ii)) / mrhi;
end

% read the parameter file
filename = strcat(dir, 'model_parameters.dat');
fid = fopen(filename, 'r');
nparam = 32;
A = fscanf(fid, '%f', [nparam+1 inf]);
fclose(fid);
A = A';
params = A(:,paramId+1);

% resampling
nresample = 2000;
resamples = zeros(nresample,1);
rands = rand(nresample,1);
for ii = 1 : nresample
   indx = find(rcdf>rands(ii), 1, 'first');
   resamples(ii) = params(indx, 4);
end

%pmin = parmin(paramId);
%pmax = parmax(paramId);
pmin = parmin(4);
pmax = parmax(4);

xii = linspace(pmin, pmax, 100001);
xi = linspace(pmin, pmax, 1001);
pdf = ksdensity(resamples, xii, 'function', 'pdf', 'support', [pmin,pmax]);
cdf = ksdensity(resamples, xi, 'function', 'cdf', 'support', [pmin,pmax]);
[~, indx] = max(pdf);
optimum = xii(indx);

ofile = ['./cdf', num2str(11), '.mat'];
save(ofile, 'xi', 'cdf');

end

