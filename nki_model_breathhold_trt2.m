function [] = nki_model_breathhold_trt2(path_folder,opt)
% generate individual model for brathhold task nki database(test retest only).
%
% SYNTAX:
% [] = HCP_IND_MODEL_BREATHHOLD_TRT2(OPT)
%
% _________________________________________________________________________
% PATH_FOLDER (string, default [pwd]) the full path to the root folder that contain the 
%   individual NKI task model variables. 
%
% OPT
%   (structure, optional) with the following fields :
%
%   TASK (string, default 'breathhold') type of trial that would be extracted, possiblr trial : 'checkerboard', 'breathhold'.
%   TRIAL_DELAY (numeric, default '1,5') time delay before each trial to estimate FIR responses.
%   TRIAL_DURATION (numeric, default '16.5') time for each trial to estimate FIR responses.
%   BASELINE_DELAY (numeric, default '2.5') time delay before each trial's baseline estimate.
%   BASELINE_DURATION (numeric, default '4') time for each trial to estimate baselie for FIR responses.
%   EXP (string, default '1400') type of TR used. Possible value '1400', '645'.

%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%

path_folder = niak_full_path (path_folder);

%% Default options
list_fields   = { 'task'       , 'trial_delay' , 'trial_duration' , 'baseline_delay' , 'baseline_duration' , 'exp'  };
list_defaults = { 'breathhold' , 1.5           ,  30.5              ,  2.5             ,  4                  , '1400' };
if (nargin > 1) && ~isempty(opt.task) 
    opt = psom_struct_defaults(opt,list_fields,list_defaults);
else
    opt = psom_struct_defaults(struct(),list_fields,list_defaults);
end

%% path and files names
data.dir_output         = path_folder;
data.name_csv_intrarun  = ['nki_model_intrarun_' lower(opt.task) '_trt2'];

%% intrarun fir
opt_csv_read.separator= sprintf(',');
data.covariates_intrarun_names = {'times','duration'};
data.covariates_intrarun_cond  = {'baseline',opt.task,'baseline',opt.task,'baseline',opt.task,'baseline',opt.task,'baseline',opt.task,'baseline',opt.task,'baseline',opt.task};
data.covariates_intrarun_values(1,1) = (144 - opt.baseline_delay) + 10;
data.covariates_intrarun_values(1,2) = opt.baseline_duration;
data.covariates_intrarun_values(2,1) = (154 - opt.trial_delay);
data.covariates_intrarun_values(2,2) = opt.trial_duration;
data.covariates_intrarun_values(3,1) = (180 - opt.baseline_delay) + 10;
data.covariates_intrarun_values(3,2) = opt.baseline_duration;
data.covariates_intrarun_values(4,1) = (190 - opt.trial_delay);
data.covariates_intrarun_values(4,2) = opt.trial_duration;
data.covariates_intrarun_values(5,1) = (216 - opt.baseline_delay) + 10;
data.covariates_intrarun_values(5,2) = opt.baseline_duration;
data.covariates_intrarun_values(6,1) = (226 - opt.trial_delay);
data.covariates_intrarun_values(6,2) = opt.trial_duration;


opt_csv_write.labels_y = data.covariates_intrarun_names;
opt_csv_write.labels_x = data.covariates_intrarun_cond;
opt_csv_write.precision = 2;
niak_write_csv(strcat(data.dir_output,data.name_csv_intrarun,'.csv'),data.covariates_intrarun_values,opt_csv_write);

