function varargout = spm_file_manager(varargin)
%__________________________________________________________________________
% Collection of functions for reading and organising data.
%
% FORMAT dat = init_dat(dir_population,S,dat)
% FORMAT modify_json_field(pth_json,field,val)
% FORMAT modify_pth_in_population(dir_population,field,npth)
%
% FORMAT help spm_file_manager>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_file_manager
    error('Not enough argument. Type ''help spm_file_manager'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id) 
    case 'init_dat'
        [varargout{1:nargout}] = init_dat(varargin{:});             
    case 'modify_json_field'
        [varargout{1:nargout}] = modify_json_field(varargin{:});             
    case 'modify_pth_in_population'
        [varargout{1:nargout}] = modify_pth_in_population(varargin{:});                 
    otherwise
        help spm_file_manager
        error('Unknown function %s. Type ''help spm_file_manager'' for help.', id)
end
%==========================================================================

%==========================================================================
function dat = init_dat(dir_population,dat)
% Reads population and meta data into a dat struct.
% FORMAT dat = init_dat(dir_population,dat)
%
% dir_population - Path to a directory containing JSON files. Each JSON
% file holds subject-specific meta data.
% S - Number of subjects to read. [S=Inf]
% dat - A cell-array that contains subject-specific data. Give as input to append 
% already exisiting dat object with new subjects or new data for already 
% existing subjects. [dat={}]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, 
    dat = {}; 
    J1  = 0; % Path to JSON file.
else
    J1  = numel(dat);
end

% Create a dictionary which will map subject names to dat indeces
dict = containers.Map;

% Get all JSON files in population directory
json_files = dir(fullfile(dir_population,'*.json'));
J          = numel(json_files);

% For displaying number of files read
base10 = floor(log10(J)) + 1;
s      = sprintf(['spm_file_manager(''init_dat'') | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],0,J);
fprintf(1, ['%-' num2str(2*base10 + 50) 's'], s);
    
tic;
for j=1:J % Loop over JSON files
    
    if ~mod(j,10)
        % Display number of files read
        fprintf(1, repmat('\b',1,2*base10 + 50));
        s = sprintf(['spm_file_manager(''init_dat'') | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],j,J);
        fprintf(1, ['%-' num2str(2*base10 + 50) 's'], s);
    end
    
    % Read subject-specific meta data
    pth_json = fullfile(dir_population,json_files(j).name);
    metadata = spm_jsonread(pth_json);       
    
    % Check and init meta_data
    metadata = check_metadata(metadata);   
    
    % Get poulation and subject name
    name       = metadata.name;
    population = metadata.population;
    
    % Create dictionary key
    key = [population '_' name];
    
    if ~dict.isKey(key)
        % Subject not in dictionary -> add subject to dictionary
        dict(key)                 = J1 + dict.Count + 1;          
        dat{dict(key)}.name       = name;
        dat{dict(key)}.population = population;
    end                    
    
    if ~isempty(metadata.modality)
        % Process imaging data
        %------------------------------------------------------------------
        modality    = metadata.modality;
        Nii         = nifti(metadata.pth);
        channel     = metadata.channel;

        if ~isfield(dat{dict(key)},'modality') || isempty(dat{dict(key)}.modality)
            % No image data exists -> create image data fields
            dat{dict(key)}.modality{1}.name = modality;
            if isempty(channel)
                % Single-channel data
                dat{dict(key)}.modality{1}.nii      = Nii;
                dat{dict(key)}.modality{1}.json.pth = pth_json;
            else
                % Multi-channel data
                dat{dict(key)}.modality{1}.channel{1}.name      = channel;
                dat{dict(key)}.modality{1}.channel{1}.nii       = Nii;
                dat{dict(key)}.modality{1}.channel{1}.json.pth  = pth_json;
            end
        else
            % Modality field already exists for subject -> append to appropriate
            % modality and channel
            M       = numel(dat{dict(key)}.modality); % Number of modalities
            has_mod = false;
            for m=1:M
                if strcmp(modality,dat{dict(key)}.modality{m}.name)
                    % Modality found -> append                               
                    if isempty(channel)
                        N = numel(dat{dict(key)}.modality{m}.nii);
                        dat{dict(key)}.modality{m}.nii(N + 1)      = Nii;
                        dat{dict(key)}.modality{m}.json(N + 1).pth = pth_json;
                    else
                        C       = numel(dat{dict(key)}.modality{m}.channel);
                        has_chn = false;
                        for c=1:C % Loop over channels
                            if strcmp(channel,dat{dict(key)}.modality{m}.channel{c}.name)
                                % Channel found -> append
                                N = numel(dat{dict(key)}.modality{m}.channel{c}.nii);
                                dat{dict(key)}.modality{m}.channel{c}.nii(N + 1)      = Nii;
                                dat{dict(key)}.modality{m}.channel{c}.json(N + 1).pth = pth_json;
                                
                                has_chn = true;
                                break
                            end                           
                        end

                        if ~has_chn
                            % Channel not found -> add channel
                            dat{dict(key)}.modality{m}.channel{C + 1}.name      = channel;
                            dat{dict(key)}.modality{m}.channel{C + 1}.nii       = Nii;
                            dat{dict(key)}.modality{m}.channel{C + 1}.json.pth  = pth_json;
                        end
                    end

                    has_mod = true;
                    break
                end
            end        

            if ~has_mod
                % Modality not found -> add modality
                if isempty(channel)
                    % Single-channel data
                    dat{dict(key)}.modality{M + 1}.nii      = Nii;
                    dat{dict(key)}.modality{M + 1}.json.pth = pth_json;
                else
                    % Multi-channel data
                    dat{dict(key)}.modality{M + 1}.channel{1}.name      = channel;
                    dat{dict(key)}.modality{M + 1}.channel{1}.nii       = Nii;
                    dat{dict(key)}.modality{M + 1}.channel{1}.json.pth  = pth_json;
                end
            end
        end
    end
    
    if ~isempty(metadata.rater)
        % Process label data
        %------------------------------------------------------------------
        rater = metadata.rater;
        Nii   = nifti(metadata.pth);
        
        if ~isfield(dat{dict(key)},'label')
            dat{dict(key)}.label{1}.name     = rater;
            dat{dict(key)}.label{1}.nii      = Nii;
            dat{dict(key)}.label{1}.json.pth = pth_json;
        else
            dat{dict(key)}.label{end + 1}.name     = rater;
            dat{dict(key)}.label{end + 1}.nii      = Nii;
            dat{dict(key)}.label{end + 1}.json.pth = pth_json;
        end
    end
    
    % Append other meta data fields (if there are any)
    %----------------------------------------------------------------------
    fn = fieldnames(metadata);
    for i=1:numel(fn)
        field_name = fn{i};
        
        if strcmp(field_name,'pth')    || ...
           strcmp(field_name,'modality') || ...
           strcmp(field_name,'name')     || ...
           strcmp(field_name,'channel')
           
            continue
        end
        
        dat{dict(key)}.(fn{i}) = metadata.(field_name);
    end            
    
    % Make sure fields are ordered alphabetically
    dat{dict(key)} = orderfields(dat{dict(key)});
end

% Display number of files read
fprintf(1, repmat('\b',1,2*base10 + 50));
s = sprintf(['spm_file_manager(''init_dat'') | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],j,J);
fprintf(1, ['%-' num2str(2*base10 + 50) 's'], s);
fprintf('\n');

fprintf('spm_file_manager(''init_dat'') | Loaded %i subjects from %s in %0.1f seconds.\n',dict.Count,dir_population,toc);
%==========================================================================

%==========================================================================
function modify_json_field(pth_json,field,val)
% Reads population and meta data into a dat struct.
% FORMAT modify_json_field(pth_json,field,val)
%
% pth_json - Path to JSON file.
% field - Field to change.
% val - Value to change to.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
a         = spm_jsonread(pth_json);
a.(field) = val;
a         = orderfields(a);
spm_jsonwrite(pth_json,a);
%==========================================================================

%==========================================================================
function modify_pth_in_population(dir_population,field,npth)
% Reads population and meta data into a dat struct.
% FORMAT modify_pth_in_population(dir_population,field,npth)
%
% dir_population - Path to a directory containing JSON files. Each JSON
% file holds subject-specific meta data.
% field - Field to change.
% npth - New path
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
json_files = dir(fullfile(dir_population,'*.json'));
J          = numel(json_files);

for s=1:J % Loop over JSON files
    pth_json = fullfile(dir_population,json_files(s).name);
        
    a           = spm_jsonread(pth_json);
    opth        = a.(field);
    ix1         = opth(1);
    [~,nam,ext] = fileparts(opth);
    if ~strcmp(ix1,filesep)
        ix1 = '';
    end
    nval        = fullfile([ix1 npth],[nam ext]);
    a.(field)   = nval;
    a           = orderfields(a);
    spm_jsonwrite(pth_json,a);
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function metadata = check_metadata(metadata)
if ~isfield(metadata,'name')
    error('~isfield(meta_data,''name'')')        
end
if ~isfield(metadata,'population')
    error('~isfield(meta_data,''population'')')        
end
if ~isfield(metadata,'modality')
    metadata.modality = '';   
end
if ~isfield(metadata,'channel')
    metadata.channel = '';   
end
if ~isfield(metadata,'rater')
    metadata.rater = '';   
end
if ~isempty(metadata.channel) && isempty(metadata.channel)
    error('~isempty(meta_data.channel) && isempty(meta_data.channel)')
end
%==========================================================================