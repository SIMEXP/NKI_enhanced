clear

%% Parameters
%stability_fir_perc_MOTORrh_hcp/stability_group/fir';
list_scale = { 'sci130_scg104_scf107' ; ...
             };
list_task = { 'rh' } ;
path_res =  '/home/pbellec/database/HCP_task/';

%% Load data
for tt = 1:length(list_task)
    path_read = [path_res 'stability_fir_perc_MOTOR' list_task{tt} '_hcp/stability_group/fir/'];
    list_files = dir([path_read 'fir_group_level_*']);
    list_files = {list_files.name};
    for ff = 1:length(list_files);
        data = load([path_read list_files{ff}],list_scale{tt});
        fir_all{tt}(:,:,ff) = data.(list_scale{tt}).fir_mean;
    end
end

%% More parameters
%list_ind = [ 19 , 33 ];
list_ind = [ 7 ];
list_color = {'r','b','g','k','p'};

%% Hierarchical clustering
hc = figure;
hp = figure;
for tt = 1:length(list_task)
    % Clustering of subtypes
    figure(hc)
    D = niak_build_distance (squeeze(fir_all{tt}(:,list_ind(tt),:)));
    hier = niak_hierarchical_clustering (-D);
    sil = niak_build_avg_silhouette(-D,hier);
    [val,nb_clust(tt)] = max(sil(1:50));
    fprintf('Task %s, max silhouette %1.2f at scale %i\n',list_task{tt},val,nb_clust(tt))
    part = niak_threshold_hierarchy (hier,struct('thresh',nb_clust(tt)));
    order = niak_hier2order (hier);
    subplot(3,length(list_task),tt)
    niak_visu_matrix(D(order,order));
    subplot(3,length(list_task),tt+length(list_task))
    niak_visu_part(part(order))
    subplot(3,length(list_task),tt+2*length(list_task))
    plot(sil)
    hold on
    plot(nb_clust(tt),val,'rx')
    
    % Show the subtypes
    figure(hp)
    subplot(1,length(list_task),tt)
    title(sprintf('Task %s',list_task{tt}));
    for cc = 1:nb_clust(tt)        
        hold on 
        plot(mean(fir_all{tt}(:,list_ind(tt),part==cc),3),list_color{cc})
    end
end