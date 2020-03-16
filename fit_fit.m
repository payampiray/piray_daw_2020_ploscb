function fit_fit(data, n, fitcat, mdls, mnames, nparams, file_hbi, hbi_tolx)
if n>length(data)
    return;
end

if nargin<7, file_hbi = ''; end
if nargin<8, hbi_tolx = 0.01; end

v0 = 6.25;

pipedir   = getdefaults('pipedir');
fitdir    = fullfile(pipedir,fitcat); makedir(fitdir);
tempdir   = getdefaults('tempdir');
lapdir    = fullfile(tempdir,fitcat,'laps'); 

N         = length(data);

%-------------------
ii = 1:length(mdls);
nn = n;

fit_hbi = 0;
if ~isempty(file_hbi) && n==0
    fit_hbi = 1;
end

%-------------------
loglik = zeros(1,length(ii));
flaps = cell(1,length(ii));

for i     = ii
    model = mdls{i};   
    mname = mnames{i};
    fname = fullfile(fitdir,sprintf('lap_%s.mat',mname));
    flaps{i} = fname;
        
    loglik(i) = model(rand(1,nparams(i)),data{1});

    if ~fit_hbi
        makedir(lapdir);
        if ~exist(fname,'file')
            d     = nparams(i);
            config.range = [-2*ones(1,d);2*ones(1,d)];
            config.numinit = 10;
            prior = struct('mean',zeros(d,1),'variance',v0);
            for n = nn
                loglik(i) = model(rand(1,d),data{n});
                flap = fullfile(lapdir,sprintf('%s_%04d.mat',mname,n));
                if exist(flap,'file') 
                    delete(flap);
                end
                if ~exist(flap,'file') && ~usejava('desktop')
                    cbm  = cbm_lap(data(n), model, prior, [], config); %#ok<NASGU>
                    save(flap,'cbm');
                end
            end
            [fnames,ok_N] = getfileordered(lapdir,sprintf('%s_%s.mat',mname,'%04d'),1:N);
            if ok_N
                cbm = cbm_lap_aggregate(fnames); %#ok<NASGU>
                save(fname,'cbm');
            end
        end
    end  
end

if fit_hbi>0    
    
    flog        = fullfile(fitdir,sprintf('hbi_%s.log',file_hbi));
    fname_hbi   = fullfile(fitdir,sprintf('hbi_%s.mat',file_hbi));
    hbiconfig   = struct('flog',flog,'tolx',hbi_tolx);
    if ~exist(fname_hbi,'file')
        cbm_hbi(data,mdls(ii),flaps,fname_hbi,hbiconfig,[]);
    end
    if exist(fname_hbi,'file')
        fname_null = fullfile(fitdir,sprintf('hbi_%s_null.mat',file_hbi));
        if ~exist(fname_null, 'file')
            cbm_hbi_null(data,fname_hbi);
        end
    end
end

end
