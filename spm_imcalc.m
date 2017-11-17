function varargout = spm_imcalc(varargin)
%__________________________________________________________________________
% Collection of tools for image calculation (gradient, suff stat, ...).
%
% FORMAT out = spm_imcalc(('name'), input)
%
% FORMAT help spm_imcalc>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_imcalc
        error('Not enough argument. Type ''help spm_imcalc'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'foo'
            [varargout{1:nargout}] = foo(varargin{:});
        otherwise
            help spm_imcalc
            error('Unknown function %s. Type ''help spm_imcalc'' for help.', id)
    end
end