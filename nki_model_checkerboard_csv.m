function [] = hcp_model_motor_csv(root_path,opt)
% generate group and individual model for motor task HCP database.
%
% SYNTAX:
% [] = HCP_MODEL_MOTOR_CSV(OPT)
%
% _________________________________________________________________________
% root_path (string, default [pwd]) the full path to the root folder that contain the 
%   HCP Preprocessed data. 
%
% OPT
%   (structure, optional) with the following fields :
%
%   TASK (string, default 'motor') type of tasks that would be extracted. Possibles tasks are: 'EMOTION',
%       'GAMBLING','LANGUAGE','MOTOR','REST','RELATIONAL','SOCIAL','WM'.
%   EXP (string, default 'hcp') type of experiment that would be extracted

%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%
%% Default path
if isempty(root_path)
  %% Setting input/output files 
  [status,cmdout] = system ('uname -n');
  server          = strtrim(cmdout);
  if strfind(server,'lg-1r') % This is guillimin
      root_path = '/gs/scratch/yassinebha/HCP/';
      fprintf ('server: %s (Guillimin) \n ',server)
      my_user_name = getenv('USER');
  elseif strfind(server,'ip05') % this is mammouth
      root_path = '/mnt/parallel_scratch_ms2_wipe_on_april_2015/pbellec/benhajal/HCP/';
      fprintf ('server: %s (Mammouth) \n',server)
      my_user_name = getenv('USER');
  else
      switch server
          case 'peuplier' % this is peuplier
          root_path = '/media/scratch2/HCP_unproc_tmp/';
          fprintf ('server: %s\n',server)
          my_user_name = getenv('USER');
          
          case 'noisetier' % this is noisetier
          root_path = '/media/database1/';
          fprintf ('server: %s\n',server)
          my_user_name = getenv('USER');
      end
  end
   
end
root_path = niak_full_path (root_path);

%% Default options
list_fields   = { 'task'    , 'exp'};
list_defaults = { 'motor  ' , 'hcp'   };
if (nargin > 1)  
    opt = psom_struct_defaults(opt,list_fields,list_defaults);
else
    opt = psom_struct_defaults(struct(),list_fields,list_defaults);
end

%% path and files names
mkdir([root_path 'fmri_preprocess_' upper(opt.task) '_' lower(opt.exp)],'EVs/');
data.dir_output         = [root_path 'fmri_preprocess_' upper(opt.task) '_' lower(opt.exp) '/EVs/'];
data.name_csv_group     = 'hcp_model_group';
data.name_csv_intrarunLR  = 'hcp_model_intrarunLR';
data.name_csv_intrarunRL  = 'hcp_model_intrarunRL';


%  %% subjects
%  data.subs =  {'HCP100307','HCP100408'};
%  
%  %% group
%  data.covariates_group_subs        = data.subs;
%  data.covariates_group_names       = {'sex','age'};
%  data.covariates_group_values(:,1) = [0;1]; % 1=f
%  data.covariates_group_values(:,2) = [30;43];
%  
%  opt_group.labels_y = data.covariates_group_names;
%  opt_group.labels_x = data.covariates_group_subs;
%  opt_group.precision = 2;
%  
%  niak_write_csv(strcat(data.dir_output,data.name_csv_group,'.csv'),data.covariates_group_values,opt_group);


%% intrarunLR fir
data.covariates_intrarunLR_names = {'times','duration'};
data.covariates_intrarunLR_cond  = {'rh','rh','baseline','baseline'};
data.covariates_intrarunLR_values(1,1) = 9.5;
data.covariates_intrarunLR_values(1,2) = 16.5; 
data.covariates_intrarunLR_values(2,1) = 129.5;
data.covariates_intrarunLR_values(2,2) = 16.5;
data.covariates_intrarunLR_values(3,1) = 8.05;
data.covariates_intrarunLR_values(3,2) = 4;
data.covariates_intrarunLR_values(4,1) = 128.5;
data.covariates_intrarunLR_values(4,2) = 4;
opt_ind.labels_y = data.covariates_intrarunLR_names;
opt_ind.labels_x = data.covariates_intrarunLR_cond;
opt_ind.precision = 2;
niak_write_csv(strcat(data.dir_output,data.name_csv_intrarunLR,'.csv'),data.covariates_intrarunLR_values,opt_ind);

%% intrarunRL fir
data.covariates_intrarunRL_names = {'times','duration'};
data.covariates_intrarunRL_cond  = {'rh','rh','baseline','baseline'};
data.covariates_intrarunRL_values(1,1) = 84.5;
data.covariates_intrarunRL_values(1,2) = 16.5; 
data.covariates_intrarunRL_values(2,1) = 160.5;
data.covariates_intrarunRL_values(2,2) = 16.5;
data.covariates_intrarunRL_values(3,1) = 83.05;
data.covariates_intrarunRL_values(3,2) = 4;
data.covariates_intrarunRL_values(4,1) = 159.5;
data.covariates_intrarunRL_values(4,2) = 4;
opt_ind.labels_y = data.covariates_intrarunRL_names;
opt_ind.labels_x = data.covariates_intrarunRL_cond;
opt_ind.precision = 2;
niak_write_csv(strcat(data.dir_output,data.name_csv_intrarunRL,'.csv'),data.covariates_intrarunRL_values,opt_ind);
