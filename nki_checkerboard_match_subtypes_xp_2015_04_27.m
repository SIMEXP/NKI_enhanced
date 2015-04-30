clear

%% Parameters
path_data = '/media/database10/nki_enhanced/';
list_run = [1400 645];
list_net = [3 4 6 7 8 10];
labels_net = { 'limb' , 'lang' , 'pcc' , 'mot' , 'mpfc' , 'att' };
nb_samps = 100;

%% Load partitions
for rr = 1:length(list_run)
    for nn = 1:length(list_net)
        path_run = [path_data 'scores_s12_run' num2str(list_run(rr)) '/output/net_' num2str(list_net(nn)) '_' labels_net{nn} '_5clusters/'];
        file_part = [path_run 'part_run_' num2str(list_run(rr)) '_net_' num2str(list_net(nn)) '.mat'];
        data = load(file_part);
        if (rr == 1) && (nn == 1)
            part = zeros(length(data.part),length(list_run),length(list_net));
        end
        part(:,rr,nn) = data.part;
    end
end

%% Measure reproducibility of subjects' partitions across test-retest
rand_net = zeros(length(list_net),1);
for nn = 1:length(list_net)
    rand_net(nn) = bmi_rand(part(:,1,nn),part(:,2,nn));
    fprintf('Rand index network %i: %1.2f\n',list_net(nn),rand_net(nn));
end

%% Test significance
samps = zeros(nb_samps,length(list_net));
for nn = 1:length(list_net)
    for num_samp = 1:nb_samps
        niak_progress(num_samp,nb_samps)
        samps(num_samp,nn) = bmi_rand(part(:,1,nn),part(randperm(size(part,1)),2,nn));
    end
end

%% Load mean s-map per subtype
[hdr,mask] = niak_read_vol([path_data 'scores_s12_run' num2str(list_run(rr)) '/mask.nii.gz']);
for rr = 1:length(list_run)
    for nn = 1:length(list_net)
        path_run = [path_data 'scores_s12_run' num2str(list_run(rr)) '/output/net_' num2str(list_net(nn)) '_' labels_net{nn} '_5clusters/'];
        file_mean = [path_run 'mean_cluster_demeaned.nii.gz'];
        [hdr,vol] = niak_read_vol(file_mean);
        if (rr == 1) && (nn == 1)
            all_mean = zeros([sum(mask(:)),size(vol,4),length(list_run),length(list_net)]);
        end
        tseries = niak_vol2tseries(vol,mask);
        all_mean(:,:,rr,nn) = tseries';
    end
end

%% Try to match test and retest
val_match = zeros(size(all_mean,2),length(list_net));
ind_match = zeros(size(all_mean,2),length(list_net));
for nn = 1:length(list_net)
    sim_subtypes = corr(all_mean(:,:,1,nn),all_mean(:,:,2,nn));
    for nn2 = 1:size(all_mean,2)
        [val,ind] = max(sim_subtypes,[],2);
        [val_max,ind_max] = max(val);
        ind = ind(ind_max);
        val_match(ind_max,nn) = val_max;
        ind_match(ind_max,nn) = ind;
        sim_subtypes(:,ind) = -Inf;
        sim_subtypes(ind_max,:) = -Inf;
    end
end