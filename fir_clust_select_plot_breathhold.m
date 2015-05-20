
function [] = fir_clust_select_plot(clust_select,path_folder,fdr_file,partition_nii_file)



%% plot selecteded nb of clusters for fir responses
%  Usage : fir_clust_select_plot(clust_select,path_folder,fdr_file,partition_nii_file)
%%%%%%%%%% example %%%%%%%%%%%%%%%%%%%%%%
%  %  clust_select=[37 12 66]                                                                           % clusters must be beteween braquet
%  %  path_folder= '/home/yassine/twins_study_basc/basc_fir/stability_fir_all_sad_blocs_test3_perc/sci140_scg140_scf147/' % the hole path folder
%  %  fdr_file='fdr_group_average_sci140_scg140_scf147'                                                  % write the file_name without .mat extension
%  %  partition_nii_file='brain_partition_threshold_group_sci140_scg140_scf147'                          % write the file_name without .nii.gz extension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
seed = psom_set_rand_seed(0);
clust_select=[ 84 , 128 , 59 , 149 , 110 , 83 ];                                                                           % clusters must be beteween braquet
path_folder= '/media/yassinebha/database2/nki_enhanced/stability_fir_perc_breathhold_1400_noscrub/stability_group/sci150_scg150_scf153_nii'; % the hole path folder
fdr_file='fdr_group_average_sci150_scg150_scf153';                                                  % write the file_name without .mat extension
partition_nii_file='brain_partition_consensus_group_sci150_scg150_scf153';                          % write the file_name without .nii.gz extension




%  read .mat file 
cd (path_folder);
load([fdr_file '.mat']);


% reorder matrix
newfir = test_fir;
newdiff = test_diff;

% fir
tmpfir.pce  = newfir.pce(:,clust_select);
tmpfir.fdr  = newfir.fdr(:,clust_select);
tmpfir.mean = newfir.mean(:,clust_select);
tmpfir.std  = newfir.std(:,clust_select);

% diff
tmpdiff.pce  = newdiff.pce(:,clust_select);
tmpdiff.fdr  = newdiff.fdr(:,clust_select);
tmpdiff.mean = newdiff.mean(:,clust_select);
tmpdiff.std  = newdiff.std(:,clust_select);

test_fir = tmpfir;
test_diff = tmpdiff;

% save matrix
str=num2str(clust_select);
str(~isstrprop(str,'alphanum')) = '_';
save( [fdr_file '_select_clust' char(str) '.mat'] , 'test_fir', 'test_diff');

%  reorder volumes
[hdr,vol] = niak_read_vol([path_folder filesep partition_nii_file '.nii.gz']);
vol2 = zeros(size(vol));
for num = 1:length(clust_select)
vol2(vol==clust_select(num)) = num;
end

% save volumes
hdr.file_name = [path_folder filesep  partition_nii_file '_select_clust' char(str) '.nii.gz'];
niak_write_vol(hdr,vol2);


% show brain patition reordered with mricron
B=size(test_fir.mean)(2);
system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_folder filesep partition_nii_file '_select_clust' char(str) '.nii.gz -m  ~/.mricron/multislice/default.ini -c NIH -l 0.02 -h ' num2str(B) ' -z & ']);

% plot fir
linewidth = 0.1;
background = [0.9 0.9 0.9]; % color of the background for non-significant responses
opt_fig.flag_legend = false;

hf = figure;
files_in = [ path_folder filesep fdr_file '_select_clust' char(str) '.mat'];
load (files_in);
[x,sizesubclust] =  size(test_fir.mean);

ymaxtmp = max(test_fir.mean) + max(test_fir.std); % added the 1.5* (avant max) because of the '*' indicating significance
ymintmp = min(test_fir.mean) - max(test_fir.std);
ymaxtmp2 = max(ymaxtmp);
ymintmp2 = min(ymintmp);
if ymaxtmp2 < 0
    ymax = 0;
else 
    ymax = ymaxtmp2 + 0.001;
end
if ymintmp2 > 0
    ymin = 0;
else 
    ymin = ymintmp2 - 0.001;
end
absmin=min((0:(size(test_fir.mean,1)+1)));
absmax=max((0:(size(test_fir.mean,1)+1)));
axisvalues = [absmin absmax];
axisvalues(1,3) = ymin;
axisvalues(1,4) = ymax;

real_absmin=min(1*(0:(size(test_fir.mean,1)-1)));
real_absmax=max(1*(0:(size(test_fir.mean,1)-1)));
xtick = [real_absmin 11 16 29 34 47 52 65 70 83 ];
ytick = [ymin ymin/2 0 ymax/4 ymax/2  ymax*3/4 ymax];
ytick = round(100*ytick)/100;


subplot(sizesubclust,sizesubclust,1+(sizesubclust-1)*sizesubclust);
ha = gca;
axis(axisvalues);
set(ha,'xtick',xtick)
set(ha,'ytick',ytick)
set(ha,'linewidth',linewidth);
for m = 1:sizesubclust
    subplot(sizesubclust,sizesubclust,m+(m-1)*sizesubclust);
    opt_fig.ind_fir = m;
    opt_fig.flag_diff = false;
    opt_fig.flag_std = true;
    opt_fig.background = background;
    opt_fig.linewidth = linewidth;
    opt_fig.axis = axisvalues;
    niak_brick_fig_fir(files_in,'',opt_fig);
    ha = gca;
    set(ha,'xtick',[xtick])
    set(ha,'ytick',[ytick])
    set(ha,'visible','on')

    for mm = 2:sizesubclust
        if mm > m
            subplot(sizesubclust,sizesubclust,mm+(m-1)*sizesubclust);
            opt_fig.flag_diff = true;
            opt_fig.flag_std = true;
            opt_fig.background = background;
            opt_fig.linewidth = linewidth;
            opt_fig.ind_fir(1,1) = m;
            opt_fig.ind_fir(1,2) = mm;
            files_out = '';
            niak_brick_fig_fir(files_in,files_out,opt_fig);
            ha = gca;
            set(ha,'xtick',xtick)
            set(ha,'ytick',ytick)
            set(ha,'visible','on')
        end
    end
end
    
    
% save .pdf and .svg files 
namesave = [ path_folder filesep 'plot_fir_select_clust' char(str)];
files_out = strcat(namesave,'.svg');
print(files_out,'-dsvg','-r600');
    
%  split cluster partition
cd ([ path_folder filesep]);
niak_brick_clusters_to_3d([partition_nii_file '_select_clust' char(str) '.nii.gz']);

%   isolated scale networks figure in mricron
for i=1:size(test_fir.mean)(2)
    if i < 10
    system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_folder filesep partition_nii_file '_select_clust' char(str) '_000' num2str(i) '.nii.gz -c jet_linear -l 0.02 -h ' num2str(B) ' -z & ']);
    system(['echo ' 'fig_' partition_nii_file '_select_clust' char(str) '_000' num2str(i)]);
    else
    system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_folder filesep partition_nii_file '_select_clust' char(str) '_00' num2str(i) '.nii.gz -c jet_linear -l 0.02 -h ' num2str(B) ' -z & ']);
    system(['echo ' 'fig_' partition_nii_file '_select_clust' char(str) '_00' num2str(i)]);
    endif
    system(['echo ' 'fig_' partition_nii_file '_select_clust' char(str) '_all']);
end


for i=1:size(test_fir.mean)(2)
    if i < 10
    system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_folder filesep partition_nii_file '_000' num2str(i) '.nii.gz -c jet_linear -l 0.02 -h ' num2str(B) ' -z & ']);
    system(['echo ' 'fig_' partition_nii_file '_s7_000' num2str(i)]);
    else
    system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_folder filesep partition_nii_file '_select_clust' char(str) '_00' num2str(i) '.nii.gz -c jet_linear -l 0.02 -h ' num2str(B) ' -z & ']);
    system(['echo ' 'fig_' partition_nii_file '_select_clust' char(str) '_00' num2str(i)]);
    endif
    
end
system(['echo ' 'fig_' partition_nii_file '_s' num2str(B) '_all']);
endfunction
