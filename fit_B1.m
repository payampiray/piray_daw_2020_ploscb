function [fname_hbi, model, fname_lap] = fit_B1(input)

bmc = 0; fitdir = 'fit_B1';
if nargin<1, input = 0; bmc=1; end

n = input(1);
% ------------

pipedir   = getdefaults('pipedir');

models    = {@fit_B_vkf, @fit_B_hgf, @fit_B_rl, @fit_B_kf};
mnames    = {'vkf','hgf','rl','kf'};
nparams   = [4 4 2 3];
file_hbi  = 'vkf_hgf_rl_kf'; hbi_tolx = 0.05;

if ~bmc
    data      = load(fullfile(pipedir,'data','dataset_matt2_noneg.mat')); data = data.data;    
    fit_fit(data, n, fitdir, models, mnames, nparams, file_hbi, hbi_tolx);
end


% model comparison
if bmc
    fname_hbi = fullfile(pipedir,fitdir,sprintf('hbi_%s.mat',file_hbi));
    M = length(mnames);
    cbm  = load(fname_hbi); 
    cbm = cbm.cbm;
    pxp  = cbm.output.protected_exceedance_prob;
    nbar = cbm.output.model_frequency;

    M0 = {'Model no', 'Model', 'PXP', 'MF'};
    M = [num2cell(1:M)' mnames' num2cell(pxp') num2cell(nbar')];
    M = [M0; M];
    
    fprintf('HBI Bayesian model comparison\n\n')    
    disp(M);
        
    fname_lap = fullfile(pipedir,fitdir,sprintf('lap_%s.mat',mnames{1}));
    model = models{1};    
    
    [fdir,filename] = fileparts(fname_lap);
    flr = fullfile(fdir,sprintf('%s_lr.mat',filename));
    lr = load(flr);
    lr = lr.mlr;    

    dlr = lr*[-1 1]';
    sdlr = nanmean(dlr>0);
    mdlr = nanmean(dlr);
    [~,plr] = ttest(dlr);
end

end
