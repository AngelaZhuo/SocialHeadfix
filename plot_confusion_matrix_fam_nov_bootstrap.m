function f = plot_confusion_matrix_fam_nov_bootstrap(d,unit_maps,ops)
%% adjusted from script Fig 1e
% MA
%% Reproduce Fig. 1e from deposited dataset
%  David Wolf, 11/2023
%
compare_classifier=1;
zscore_data =1;
% build spike count population vector
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);
time_base = -ops.pre:ops.binsize:ops.post;
time_base(end) = [];
response_bins = find(round(time_base,3)==ops.response(1)):find(round(time_base,3)==ops.response(2))-1;

% the data.(odor) contains the binned spike counts (100 ms bins) for the
% different odors and trials (dimensions: unit x time x trials)
num_its = 100;
num_units = 75;

shuffle_labels = false;

% confusion_matrix = NaN(size(pop_vec,2),size(pop_vec,2),num_its);
% tmp = cellfun(@(x) squeeze(mean(x(:,response_bins,:),2)), pop_vec, 'UniformOutput', false);
% X = cat(2,tmp{1,:})';
% y = repelem(1:size(pop_vec,2), size(pop_vec{1},3))';
% % ### there is a bug resulting in non-sense results for AON CTRL
% % bootstrapStats = bootstrapClassifierValidation(X, y); 
% % bootstrapStats = bootstrapClassifierLeaveOneOut(X, y);
% % ###

% original 
 [confusion_matrix(:,:,1), correct_or_not(:,1)] = odor_response_classifier(pop_vec,shuffle_labels,zscore_data,ops);

% bootstrap
parfor nx = 1:num_its
    
    bootIndices = randi(numel(unit_maps), numel(unit_maps), 1);
%     rand_units = randperm(numel(unit_maps));
%     curr_units = (rand_units(1:num_units));
    curr_pop_vec = cellfun(@(x) x(bootIndices,:,:), pop_vec, 'UniformOutput', false);
    [confusion_matrix(:,:,nx+1), correct_or_not(:,nx+1)] = odor_response_classifier(curr_pop_vec,shuffle_labels,zscore_data,ops);
end
%% compare to classifier trained with shuffled labels
if compare_classifier
    num_it = 10;
   for i = 1:num_it
       shuffle_labels = true;
       [~, correct_or_not_shuffle] = odor_response_classifier(pop_vec,shuffle_labels,zscore_data,ops);
       contingency_table = table([nnz(correct_or_not(:,1));nnz(correct_or_not_shuffle)],...
           [nnz(~correct_or_not(:,1));nnz(~correct_or_not_shuffle)],'VariableNames',{'correct','not'},'RowNames',{'true','shuffle'});
       [h,p(i),stats] = fishertest(contingency_table);
   end
end

%% plot
f = figure('Position',[1 1 30 10]); tiledlayout(1,3);
if size(confusion_matrix,1)>3
    odor_names = ["ylang ylang"  "peanut butter" "BL6 #1" "BL6 #2" "CD1"];
else
    odor_names = ["fam" "nov1" "nov2"];
end
curr_odors = odor_names(1:size(confusion_matrix,1));
%%
for hx = 1:2
nexttile
if hx ==1
    heatmap(curr_odors,curr_odors,(confusion_matrix(:,:,1)),'GridVisible','off');
    accuracy = mean(correct_or_not(:,1));
else
    heatmap(curr_odors,curr_odors,median(confusion_matrix(:,:,2:end),3),'GridVisible','off');
    accuracy = mean(mean(correct_or_not(:,2:end),1));
end
% imagesc(confusion_matrix);
ylabel('actual');
xlabel('predicted');
if compare_classifier && hx == 1
    title({['accuracy ' num2str(accuracy)];[' p:' num2str(median(p)) ' tested against ' num2str(num_it) ' shuffles']});
else
    title({'median across bootstraps' ;['accuracy ' num2str(accuracy)]});
end
drawnow
colorbar
axis = gca;
if contains(class(axis),'heatmap','IgnoreCase',true)
    axis.ColorLimits= [0 size(pop_vec{1, 1},3)];
else
    axis.XTick = 1:size(confusion_matrix,1);
    axis.YTick = 1:size(confusion_matrix,1);
    axis.XTickLabel = odor_names(1:size(confusion_matrix,1));
    axis.YTickLabel = odor_names(1:size(confusion_matrix,1));
    set(get(axis, 'XLabel'), 'FontSize', 6);
    set(get(axis, 'YLabel'), 'FontSize', 6);
    set(axis, 'FontSize', 6);
    set(axis, 'FontName', 'Arial');
    set(get(axis, 'Title'), 'FontSize', 8);
    axis.CLim = [0 size(pop_vec{1, 1},3)];
end

colormap('pink');
end

nexttile
histogram(mean(correct_or_not,1))
title(['accuracy distribution ' num2str(num_its)  ' iterations'])

% f.Units = 'centimeters';
% f.Position = [3 3 8.3 5.6];




