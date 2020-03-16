function [filelist,ok] = getfileordered(fdir,filter,order)

generr = 1;
if nargout>1, generr = 0; end

% filter such as 'job_temp%02d*.mat' or 'job_temp%02d.mat'

currdir = pwd;

filelist = cell(length(order),1);
oks = nan(length(order),1);
cd(fdir);
for i=1:length(order)
    oi = order(i);        
    fi = sprintf(filter,oi);
    fname = dir(fi);
    
    oks(i) = length(fname);
    if generr
        if oks(i)~=1
            error('filter does not correspond to 1 and only 1 file!');
        end
    end
    
    if oks(i)==1
        filelist{i} = fullfile(fdir,fname.name);
    end
    
end
cd(currdir);

ok = all(oks==1);

end