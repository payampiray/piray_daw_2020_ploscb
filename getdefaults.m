function varargout = getdefaults(action,varargin)

st = dbstack('-completenames');
maindir = (fileparts(st(1).file));
pipedir = fullfile(maindir,'pipe');
tempdir = fullfile(maindir,'temp');

switch action
case 'pipedir'
    varargout{:} = pipedir;
    addpath(fullfile(maindir,'tools'));
    addpath(fullfile(maindir,'cbm_local'));    
    
case 'tempdir'
    varargout{:} = tempdir;
    
% figure properties
    case 'alf'
        varargout{:} = .6;
    case 'fs'
        varargout{:} = 12;
    case 'fsy'
        varargout{:} = 14;
    case 'fsl'
        varargout{:} = 14;
    case 'fst'
        varargout{:} = 18;
    case 'fsA'
        varargout{:} = 18;
    case 'fn'
        varargout{:} = 'Arial';
    case 'xsA'        
        varargout{:} = -.15;
    case 'ysA'
        varargout{:} = 1.1;
    case 'fnA'
        varargout{:} = 'Arial';
        
        
    case 'yst'
        varargout{:} = 1.15;        
    case 'fnt'
        varargout{:} = 'Helvetica';
        varargout{:} = 'Arial';
        
    case 'abc'
        varargout{:} = 'ABCDEF';
    case 'colmap'
        varargout{:} = [1 .4 0;0 .4 1];
    case 'cmaphbi'
        varargout{:} = [1 .2 .2];
        
otherwise
error('%s does not exist',action);
end
end

