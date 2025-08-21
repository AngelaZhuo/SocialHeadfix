function euclidian = compute_popVec_distance_from_baseline(pop_vec,ops)
%% Computes Population Vector for response window between conditions
%
% Response window population vector: counts spikes in response window
%
% Output:
%   euclidian distance to baseline 
%%

% prepare timebase
% time_base = round(-ops.pre:ops.binsize:ops.post,-log10(ops.binsize));
time_base = round(-ops.pre:ops.binsize:ops.post,3);
time_base(end) = [];
%
% % get baseline bins
baseline_bins = find(ops.baseline(1)==time_base):find(ops.baseline(2)==time_base)-1;

% buid baseline vector
base_vector = cell(1,size(pop_vec,2));
cum_base_vec = [];
for ii = 1:size(pop_vec,2)
    % baseline_vector{1,ii} = nan(size(pop_vec{1,ii},[1 3]));
    % for tri = 1:size(pop_vec{1,ii},3)
    %     baseline_vector{1,ii}(:,tri) = mean(pop_vec{1,ii}(:,baseline_bins,tri),2); %units x trials
    % end
    base_vector{1,ii} = squeeze(mean(pop_vec{1,ii}(:,baseline_bins,:),2)); 
    cum_base_vec = [cum_base_vec base_vector{1,ii}]; %units x trials
end
average_base_vec = mean (cum_base_vec,2);

euclidian = cell(2,size(pop_vec,2));
for ii = 1:size(euclidian,2) %conditions
    euclidian{1,ii} = nan(size(pop_vec{1,ii},3),size(pop_vec{1,ii},2));
    euclidian{2,ii} = nan(size(pop_vec{1,ii},3),size(pop_vec{1,ii},2));
    for tri = 1:size(pop_vec{1,ii},3) %trials
        euclidian{1,ii}(tri,:) = pdist2(average_base_vec',pop_vec{1,ii}(:,:,tri)')./sqrt(size(pop_vec{1,ii},1));
        euclidian{2,ii}(tri,:) = pdist2(base_vector{1,ii}(:,tri)',pop_vec{1,ii}(:,:,tri)')./sqrt(size(pop_vec{1,ii},1));
%         euclidian{ii}(tri,:) = pdist2(average_base_vec',pop_vec{1,ii}(:,:,tri)');
    end
end

end