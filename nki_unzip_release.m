%% script to unzip compressed raw release from nki enhanced database

%set pth in and out
root_path = '/media/database4/nki_enhanced/';
path_out  = '/media/scratch2/NKI_enhanced/raw_mnc_all/';
% grab a specofic set of release
nb_release = 0;
list_release_raw = dir([root_path 'release*' ]);
for num_r = 1:length(list_release_raw)
    if ~ismember(list_release_raw(num_r).name,{'.','..','octave-workspace','octave_core','release1','release2'})
       nb_release = nb_release + 1;
       sprintf('Adding %s', list_release_raw(num_r).name)
       list_release{nb_release} = list_release_raw(num_r).name;     
    else 
       sprintf('%s is discarded', list_release_raw(num_r).name)
    end  
end   
% uncompresse them in the output folder
for num_a = 1:length(list_release)
    path_release = [ root_path list_release{num_a} filesep 'raw_nii_compressed/' ];
    path_tmp = [ root_path list_release{num_a} ]
    list_archive = dir ([ path_release 'group_*']);
    for num_g = 1:length(list_archive)
        command1 = sprintf('tar -xvf %s%s  -C %s/. ',path_release,list_archive(num_g).name,path_tmp);
        fprintf('Unpacking archive number %i from %i \n',num_g,length(list_archive))
        system(command1);
        cammand2 = sprintf('scp -rv  %s/group_*/* %s.',path_tmp,path_out)
    end
end