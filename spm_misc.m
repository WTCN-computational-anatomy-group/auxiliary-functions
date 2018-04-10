function varargout = spm_misc(varargin)
%__________________________________________________________________________
% Collection of miscellaneous functions.
%
% FORMAT manage_parpool(num_workers)
% FORMAT nw = nbr_parfor_workers
%
% FORMAT help spm_parfor>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_parfor
        error('Not enough argument. Type ''help spm_parfor'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'manage_parpool'
            [varargout{1:nargout}] = manage_parpool(varargin{:});
        case 'nbr_parfor_workers'
            [varargout{1:nargout}] = nbr_parfor_workers(varargin{:});            
        otherwise
            help spm_parfor
            error('Unknown function %s. Type ''help spm_parfor'' for help.', id)
    end
end
%==========================================================================

%==========================================================================
function manage_parpool(num_workers)
% Start/stop parallel pool
% FORMAT manage_parpool(num_workers)
% num_workers - Number of parfor workers
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
    nw = spm_misc('nbr_parfor_workers');
    if num_workers>nw
        num_workers = nw;
    end

    poolobj = gcp('nocreate');
    if ~isempty(poolobj) && num_workers==0
        delete(poolobj);
    elseif ~isempty(poolobj) && poolobj.NumWorkers~=num_workers
        delete(poolobj);
        parpool('local',num_workers);
    elseif isempty(poolobj) && num_workers
        parpool('local',num_workers);
    end
end
%==========================================================================

%==========================================================================
function nw = nbr_parfor_workers
% Get number of CPU cores
% FORMAT nw = nbr_parfor_workers
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
    c  = parcluster('local');
    nw = c.NumWorkers;
end
%==========================================================================
