% copy unprocessed HCP task data from connectome in box hard drives to a remote server 
 
path_data = '/media/scratch2/NKI_enhanced/raw_mnc_all/';
path_out  = '/sb/project/gsf-624-aa/database/NKI_enhanced/raw_task/';
server    = 'guillimin';
user_name = 'yassinebha';

% Grab subjects list
list_subject_raw = dir(path_data);
nb_subject = 0;
for num_ss = 1:length(list_subject_raw)
    if ~ismember(list_subject_raw(num_ss).name,{'.','..'}) && isdir([ path_data list_subject_raw(num_ss).name filesep 'TfMRI_breathHold_1400/' ])...
    && isdir([ path_data list_subject_raw(num_ss).name filesep 'TfMRI_visualCheckerboard_1400/' ]) && isdir([ path_data list_subject_raw(num_ss).name filesep 'TfMRI_visualCheckerboard_645/' ])
      nb_subject = nb_subject + 1;
      sprintf('Adding subject %s', list_subject_raw(num_ss).name)
      list_subject{nb_subject} = list_subject_raw(num_ss).name;     
    else 
      sprintf('subject %s is discarded', list_subject_raw(num_ss).name)
    end  
end

% Rsync folders
for nn=1:length(list_subject);
sprintf('Syncing subjet %s ', list_subject{nn})
system(['ssh ' user_name '@' server ' mkdir -p ' path_out list_subject{nn} filesep]);
system(['rsync -ravv -f"- session_1/" ' path_data list_subject{nn} filesep ' ' user_name '@' server ':' path_out list_subject{nn} filesep]);   
end
