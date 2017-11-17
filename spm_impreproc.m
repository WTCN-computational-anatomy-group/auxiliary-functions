function varargout = spm_impreproc(varargin)
%__________________________________________________________________________
% Collection of tools for image preprocessing.
%
% FORMAT out = spm_impreproc(('name'), input)
%
% FORMAT help spm_impreproc>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_impreproc
        error('Not enough argument. Type ''help spm_impreproc'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'foo'
            [varargout{1:nargout}] = foo(varargin{:});
        otherwise
            help spm_impreproc
            error('Unknown function %s. Type ''help spm_impreproc'' for help.', id)
    end
end
