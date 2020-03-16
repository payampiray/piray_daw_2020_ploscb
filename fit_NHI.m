function M = fit_NHI(lme,mnames)
K = length(mnames);
[~,nbar0,~,pxp0] = spm_BMS(lme);
% F = sum(lme); F = F-max(F);

M = [num2cell(1:K)' mnames' num2cell(pxp0') num2cell(nbar0')];

M0 = {'Model no', 'Model', 'PXP', 'MF'};
M = [M0; M];
fprintf('Bayesian model comparison\n\n')
disp(M);
end