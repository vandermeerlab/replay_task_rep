function set_replay_path(varargin)
% function set_replay_path(varargin)
%
% sets path for the project, and returns location of
% data folders to analyze
%
% WARNING: restores default path
%
% INPUTS
%
% (optional cfg struct fields, values below are defaults)\
%   cfg.vandermeerlab = 1; % add vandermeerlab codebase
%
% OUTPUTS
%
% (none)

cfg.vandermeerlab = 1;

if nargin == 1
    cfg = varargin{1};
elseif nargin > 1
    error('set_replay_path() requires 0 or 1 input arguments.');
end
restoredefaultpath;

if ispc
    machinename = getenv('COMPUTERNAME');
    filesep = '\';
elseif ismac
    machinename = getenv('USER');
    filesep = '/';
else
    machinename = getenv('HOSTNAME');
    filesep = '/';
end

% get base file path where repo lives
switch machinename
    case {'mac'} % add case for your machine
        base_fp = '/Users/mac/Projects/replay_task_rep/analysis';
    case {'PROMETHEUS'}
        base_fp = 'C:\Users\mvdmlab\Documents\GitHub\replay_task_rep\analysis';
end

if cfg.vandermeerlab
    addpath(genpath(cat(2,'..',filesep,'vandermeerlab',filesep,'code-matlab',filesep,'shared')));
    addpath(genpath(cat(2,'..',filesep,'vandermeerlab',filesep,'code-matlab',filesep,'tasks',filesep,'Alyssa_Tmaze')));
    addpath(genpath(cat(2,'..',filesep,'vandermeerlab',filesep,'code-matlab',filesep,'graveyard')));
end

% add to path
addpath(genpath(cat(2,base_fp,filesep,'utils')));