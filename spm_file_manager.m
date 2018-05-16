function varargout = spm_file_manager(varargin)
%__________________________________________________________________________
% Collection of functions for reading and organising data.
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
    otherwise
        help spm_file_manager
        error('Unknown function %s. Type ''help spm_file_manager'' for help.', id)
end
%==========================================================================

%==========================================================================
function dat = init_dat(dir_population,S,dat)
% Reads population and meta data into a dat struct.
% FORMAT dat = init_dat(dir_population,S,dat)
%
% dir_population - Path to a directory containing JSON files. Each JSON
% file holds subject-specific meta data.
% S - Number of subjects to read. [S=Inf]
% dat - A struct that contains subject-specific data. Give as input to append 
% already exisiting dat structure with new subjects or new data for already 
% existing subjects. [dat=struct]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, S    = Inf; end
if nargin<3, 
    dat = struct; 
    S1  = 0;
else
    S1  = numel(dat);
end

% Create a dictionary which will map subject names to dat indeces
dict = containers.Map;

% Get all JSON files in population directory
json_files = dir(fullfile(dir_population,'*.json'));

S0 = numel(json_files);
S  = min(S0,S);
for s=1:S % Loop over subjects in population
    
    % Read subject-specific meta data
    pth_json  = fullfile(dir_population,json_files(s).name);
    meta_data = spm_jsonread(pth_json);       
    
    % Check and init meta_data
    meta_data = check_meta_data(meta_data);   
    
    % Get subject name
    name = meta_data.name;
    
    if ~dict.isKey(name)
        % Subject not in dictionary -> add subject to dictionary
        dict(name)           = S1 + dict.Count + 1;          
        dat(dict(name)).name = name;
    end                    
    
    if ~isempty(meta_data.modality)
        % Get imaging data
        %------------------------------------------------------------------
        modality    = meta_data.modality;
        Nii         = nifti(meta_data.pth);
        channel     = meta_data.channel;

        if ~isfield(dat(dict(name)),'modality') || isempty(dat(dict(name)).modality)
            % No image data exists -> create image data fields
            dat(dict(name)).modality.name = modality;
            if isempty(channel)
                % Single-channel data
                dat(dict(name)).modality.nii      = Nii;
                dat(dict(name)).modality.json.pth = pth_json;
            else
                % Multi-channel data
                dat(dict(name)).modality.channel.name      = channel;
                dat(dict(name)).modality.channel.nii       = Nii;
                dat(dict(name)).modality.channel.json.pth  = pth_json;
            end
        else
            % Modality field already exists for subject -> append to appropriate
            % modality and channel
            M       = numel(dat(dict(name)).modality); % Number of modalities
            has_mod = false;
            for m=1:M
                if strcmp(modality,dat(dict(name)).modality(m).name)
                    % Modality found -> append                               
                    if isempty(channel)
                        N = numel(dat(dict(name)).modality(m).nii);
                        dat(dict(name)).modality(m).nii(N + 1)      = Nii;
                        dat(dict(name)).modality(m).json(N + 1).pth = pth_json;
                    else
                        C       = numel(dat(dict(name)).modality(m).channel);
                        has_chn = false;
                        for c=1:C % Loop over channels
                            if strcmp(channel,dat(dict(name)).modality(m).channel(c).name)
                                % Channel found -> append
                                N = numel(dat(dict(name)).modality(m).channel(c).nii);
                                dat(dict(name)).modality(m).channel(c).nii(N + 1)      = Nii;
                                dat(dict(name)).modality(m).channel(c).json(N + 1).pth = pth_json;
                            end

                            has_chn = true;
                            break
                        end

                        if ~has_chn
                            % Channel not found -> add channel
                            dat(dict(name)).modality(m).channel(C + 1).name      = channel;
                            dat(dict(name)).modality(m).channel(C + 1).nii       = Nii;
                            dat(dict(name)).modality(m).channel(C + 1).json.pth  = pth_json;
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
                    dat(dict(name)).modality(M + 1).nii      = Nii;
                    dat(dict(name)).modality(M + 1).json.pth = pth_json;
                else
                    % Multi-channel data
                    dat(dict(name)).modality(M + 1).channel.name      = channel;
                    dat(dict(name)).modality(M + 1).channel.nii       = Nii;
                    dat(dict(name)).modality(M + 1).channel.json.pth  = pth_json;
                end
            end
        end
    end
    
    % Append other meta data fields (if there are any)
    %----------------------------------------------------------------------
    fn = fieldnames(meta_data);
    for i=1:numel(fn)
        field_name = fn{i};
        
        if strcmp(field_name,'pth')    || ...
           strcmp(field_name,'modality') || ...
           strcmp(field_name,'name')     || ...
           strcmp(field_name,'channel')
           
            continue
        end
        
        dat(dict(name)).(fn{i}) = meta_data.(field_name);
    end            
end

fprintf('Loaded %i subjects into dat from %s.\n',S,dir_population);
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function meta_data = check_meta_data(meta_data)
if ~isfield(meta_data,'name')
    error('~isfield(meta_data,''name'')')        
end
if ~isfield(meta_data,'modality')
    meta_data.modality = '';   
end
if ~isfield(meta_data,'channel')
    meta_data.channel = '';   
end
if ~isempty(meta_data.channel) && isempty(meta_data.channel)
    error('~isempty(meta_data.channel) && isempty(meta_data.channel)')
end
%==========================================================================