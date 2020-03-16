function fit_B2(input)
% Comparing two version of hgf (in hgf2, omega if fixed at -4)

bmc = 0; fitdir = 'fit_B1';
if nargin<1, input = 0; bmc=1; end

n = input(1);
run(n,fitdir,bmc);
end

function run(n,fitdir,bmc)
pipedir   = getdefaults('pipedir');

models    = {@fit_B_hgf,@fit_B_hgf2};
mnames    = {'hgf','hgf2'};
nparams   = [4 3];

if ~ bmc
    data      = load(fullfile(pipedir,'data','dataset_matt2_noneg.mat')); data = data.data;    
    fit_fit(data, n, fitdir, models, mnames, nparams);
end

if bmc
    M =length(mnames);
    lme = cell(1,M);
    for i = 1:length(mnames)
        fname = fullfile(pipedir, fitdir, sprintf('lap_%s.mat', mnames{i}) );
        cbm   = load(fname); 
        cbm = cbm.cbm;
        lme{i}  = cbm.output.log_evidence;
    end   
    lme = cell2mat(lme);
    fit_NHI(lme,mnames);
end
end