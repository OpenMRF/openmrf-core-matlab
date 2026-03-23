function fullpath = save_study_results(res, varargin)
%SAVE_STUDY_RESULTS Save reconstruction outputs to a MAT file.
%
% Usage:
%   save_study_results(res)
%   save_study_results(res, params)
%   save_study_results(res, params, study_info)
%   save_study_results(res, params, study_info, study_name, default_dir, use_file_dialog)
%   save_study_results(res, 'use_file_dialog', false)
%   save_study_results(res, params, 'use_file_dialog', false)
%
% Name-value options:
%   'params'
%   'study_info'
%   'study_name'
%   'default_dir'
%   'use_file_dialog' (or 'prompt')
%
% use_file_dialog / prompt:
%   true  -> open file explorer dialog (default)
%   false -> save directly if both directory and file name are available

    if isempty(res)
        error('save_study_results:MissingRes', 'Input ''res'' must not be empty.');
    end

    opts = struct();
    opts.params = struct();
    opts.study_info = struct();
    opts.study_name = '';
    opts.default_dir = '';
    opts.use_file_dialog = true;

    [opts, arg_idx] = parse_positional_inputs(opts, varargin);
    opts = parse_name_value_inputs(opts, varargin(arg_idx:end));

    params = opts.params;
    study_info = opts.study_info;
    study_name = opts.study_name;
    default_dir = opts.default_dir;
    use_file_dialog = opts.use_file_dialog;

    if isempty(study_name)
        study_name = get_study_name(study_info);
    end
    if isempty(default_dir)
        default_dir = get_study_path(study_info);
    end

    savename = '';
    if ~isempty(study_name)
        [~, base] = fileparts(char(study_name));
        if ~isempty(base)
            savename = [base '.mat'];
        end
    end

    if use_file_dialog
        if ~isempty(default_dir) && ~isempty(savename)
            default_path = fullfile(default_dir, savename);
        elseif ~isempty(default_dir)
            default_path = default_dir;
        elseif ~isempty(savename)
            default_path = savename;
        else
            default_path = '';
        end

        [filename, pathname] = uiputfile('*.mat', 'Save results as', default_path);
        if isequal(filename, 0) || isequal(pathname, 0)
            disp('Save canceled by user.');
            fullpath = '';
            return;
        end
        fullpath = fullfile(pathname, filename);
    else
        if isempty(default_dir) || isempty(savename)
            disp('Nothing saved: output path or file name could not be determined.');
            fullpath = '';
            return;
        end
        if ~exist(default_dir, 'dir')
            disp('Nothing saved: output directory does not exist.');
            fullpath = '';
            return;
        end
        fullpath = fullfile(default_dir, savename);
    end

    save(fullpath, 'res', 'params', 'study_info', '-v7.3');
    fprintf('Saved results to: %s\n', fullpath);
end

function [opts, arg_idx] = parse_positional_inputs(opts, args)
    arg_idx = 1;

    if arg_idx <= numel(args) && ~is_option_name(args{arg_idx})
        opts.params = args{arg_idx};
        arg_idx = arg_idx + 1;
    end
    if arg_idx <= numel(args) && ~is_option_name(args{arg_idx})
        opts.study_info = args{arg_idx};
        arg_idx = arg_idx + 1;
    end
    if arg_idx <= numel(args) && ~is_option_name(args{arg_idx})
        opts.study_name = char(args{arg_idx});
        arg_idx = arg_idx + 1;
    end
    if arg_idx <= numel(args) && ~is_option_name(args{arg_idx})
        opts.default_dir = char(args{arg_idx});
        arg_idx = arg_idx + 1;
    end
    if arg_idx <= numel(args) && ~is_option_name(args{arg_idx})
        opts.use_file_dialog = logical(args{arg_idx});
        arg_idx = arg_idx + 1;
    end
end

function opts = parse_name_value_inputs(opts, args)
    if isempty(args)
        return;
    end
    if mod(numel(args), 2) ~= 0
        error('save_study_results:InvalidInput', 'Name-value inputs must come in pairs.');
    end

    for k = 1:2:numel(args)
        key = lower(char(args{k}));
        value = args{k + 1};
        switch key
            case 'params'
                opts.params = value;
            case {'study_info', 'studyinfo'}
                opts.study_info = value;
            case {'study_name', 'studyname'}
                opts.study_name = char(value);
            case {'default_dir', 'defaultdir'}
                opts.default_dir = char(value);
            case {'use_file_dialog', 'usefiledialog', 'prompt'}
                opts.use_file_dialog = logical(value);
            otherwise
                error('save_study_results:InvalidOption', 'Unknown option: %s', key);
        end
    end
end

function tf = is_option_name(value)
    tf = false;
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        return;
    end
    key = lower(char(value));
    tf = any(strcmp(key, {'params', 'study_info', 'studyinfo', 'study_name', 'studyname', ...
                          'default_dir', 'defaultdir', 'use_file_dialog', 'usefiledialog', 'prompt'}));
end

function study_name = get_study_name(study_info)
    study_name = '';
    if isstruct(study_info) && isfield(study_info, 'study_name') && ~isempty(study_info.study_name)
        study_name = char(study_info.study_name);
    end
end

function study_path = get_study_path(study_info)
    study_path = '';
    if isstruct(study_info) && isfield(study_info, 'study_path') && ~isempty(study_info.study_path)
        study_path = char(study_info.study_path);
    end
end
