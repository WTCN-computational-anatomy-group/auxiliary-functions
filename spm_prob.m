function varargout = spm_prob(varargin)
%__________________________________________________________________________
% Collection of tools for probability distributions (PDF, KL-div, ...).
%
% FORMAT out = spm_prob(('name'), input)
%
% FORMAT help spm_prob>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_prob
        error('Not enough argument. Type ''help spm_prob'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'foo'
            [varargout{1:nargout}] = foo(varargin{:});
        otherwise
            help spm_prob
            error('Unknown function %s. Type ''help spm_prob'' for help.', id)
    end
end

