function [euclidian, cosine] = compute_popVec_distance_between_condition(pop_vec,ops)
%% Computes Population Vector for response window between conditions
%
% Response window population vector: counts spikes in response window
%
% Output:
%   euclidian/cosine distance between conditions. Upper triangle spike
%   count, lower triangle normalized spike count
%%

% euclidian_to_base = compute_popVec_distance_from_baseline(pop_vec,ops);

euclidian = repmat({nan(size(pop_vec{1,1},3),size(pop_vec{1,1},3),size(ops.response,1))},size(pop_vec,2),size(pop_vec,2));
cosine = euclidian;

for rx = 1:size(ops.response,1)
% prepare timebase
time_base = round(-ops.pre:ops.binsize:ops.post,3);
time_base(end) = [];
% 

% % get baseline bins
baseline_bins = find(ops.baseline(1)==time_base):find(ops.baseline(2)==time_base)-1;

% get response window bins
response_bins = find(ops.response(rx,1)==time_base):find(ops.response(rx,2)==time_base)-1;


% buid response vector 
resp_vector = cell(2,size(pop_vec,2));
baseline_vector = cell(1,size(pop_vec,2));
for ii = 1:size(pop_vec,2)
   for tr = 1:size(pop_vec{1,ii},3)
%        resp_vector{1,ii}(:,tr) = mean(pop_vec{1,ii}(:,response_bins,tr),2); %units x trials
%        baseline_vector{1,ii}(:,tr) = mean(pop_vec{1,ii}(:,baseline_bins,tr),2); %units x trials
       resp_vector{1,ii}(:,tr) = sum(pop_vec{1,ii}(:,response_bins,tr),2); %units x trials
       baseline_vector{1,ii}(:,tr) = sum(pop_vec{1,ii}(:,baseline_bins,tr),2); %units x trials
   end
end

for ii = 1:size(pop_vec,2)
%     for tr = 1:size(pop_vec{1,ii},3)
%         resp_vector{2,ii}(:,tr) = sum(pop_vec{2,ii}(:,response_bins,tr),2); %units x trials
%     end
   resp_vector{2,ii} = (resp_vector{1,ii}-mean(baseline_vector{1,ii},2))./std(baseline_vector{1,ii},0,2); % zscore
%    resp_vector{2,ii} = (resp_vector{1,ii}-mean([baseline_vector{1,:}],2))./std([baseline_vector{1,:}],0,2);
% % % normalize response vector by distance to baseline
%     for tr = 1:size(pop_vec{1,ii},3)
%         resp_vector{2,ii}(:,tr) = mean(pop_vec{1,ii}(:,response_bins,tr)./euclidian_to_base{1,ii}(tr,response_bins),2); %units x trials
%     end
end
nan_log = any(isnan([resp_vector{2,:}]),2);
for ii = 1:size(pop_vec,2)
    resp_vector{2,ii}(nan_log,:) = [];
end



% distance between trials
for ii = 1:size(euclidian,1) %conditions
    for jj = 1:size(euclidian,2) %conditions
%         euclidian{ii,jj} = nan(size(pop_vec{1,ii},3),size(pop_vec{1,ii},3));
%         cosine{ii,jj} = nan(size(pop_vec{1,ii},3),size(pop_vec{1,ii},3));
        for tri = 1:size(pop_vec{1,ii},3) %trials
            for trj = 1:size(pop_vec{1,ii},3) %trials

                % euclidian distance: lower triangle is normalized                
                if ii<jj
                    euclidian{ii,jj}(tri,trj,rx) = pdist2(resp_vector{1,ii}(:,tri)',resp_vector{1,jj}(:,trj)');%/sqrt(size(resp_vector{1,ii},1));
                elseif ii==jj %diagonal = within-odor distance
                    if tri<trj
                        euclidian{ii,jj}(tri,trj,rx) = pdist2(resp_vector{1,ii}(:,tri)',resp_vector{1,jj}(:,trj)');%/sqrt(size(resp_vector{1,ii},1));
                    elseif tri>trj
                        euclidian{ii,jj}(tri,trj,rx) = pdist2(resp_vector{2,ii}(:,tri)',resp_vector{2,jj}(:,trj)')./sqrt(size(resp_vector{2,ii},1));
                    end
                else
                    euclidian{ii,jj}(tri,trj,rx) = pdist2(resp_vector{2,ii}(:,tri)',resp_vector{2,jj}(:,trj)')./sqrt(size(resp_vector{2,ii},1));
                end
                % cosine distance
                if ii<jj
                    cosine{ii,jj}(tri,trj,rx) = pdist2(resp_vector{1,ii}(:,tri)',resp_vector{1,jj}(:,trj)','cosine');
                elseif ii==jj %diagonal = within-odor distance
                    if tri<trj
                        cosine{ii,jj}(tri,trj,rx) = pdist2(resp_vector{1,ii}(:,tri)',resp_vector{1,jj}(:,trj)','cosine');
                    elseif tri>trj
                        cosine{ii,jj}(tri,trj,rx) = pdist2(resp_vector{2,ii}(:,tri)',resp_vector{2,jj}(:,trj)','cosine');
                    end
                else
                    cosine{ii,jj}(tri,trj,rx) = pdist2(resp_vector{2,ii}(:,tri)',resp_vector{2,jj}(:,trj)','cosine');
                end
            end
        end
    end
end
end
euclidian = cellfun(@(x) mean(x, 3), euclidian, 'UniformOutput', false);
cosine = cellfun(@(x) mean(x, 3), cosine, 'UniformOutput', false);
end