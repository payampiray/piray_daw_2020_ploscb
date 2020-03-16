function [fname_hbi, model, fname_lap] = fit_A1(input)

bmc = 0; fitdir = 'fit_A1';
if nargin<1, input = 0; bmc=1; end

n = input(1);
% ------------


pipedir   = getdefaults('pipedir');

models    = {@fit_A_vkf, @fit_A_hgf, @fit_A_rl, @fit_A_kf};
mnames    = {'vkf','hgf','rl','kf'};
nparams   = [7 7 5 6];
file_hbi  = 'vkf_hgf_rl_kf';

if ~bmc
    data      = load(fullfile(pipedir,'data','dataset1.mat')); data = data.data;    
    fit_fit(data, n, fitdir, models, mnames, nparams, file_hbi);
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
    
    fprintf('HBI Bayesian model comparison\n\n');
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
