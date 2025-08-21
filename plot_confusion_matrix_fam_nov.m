function [f,p] = plot_confusion_matrix_fam_nov(d,unit_maps,ops)
%% adjusted from script Fig 1e
% MA
%% Reproduce Fig. 1e from deposited dataset
%  David Wolf, 11/2023
%
compare_classifier=1;
zscore_data =1;
% build spike count population vector
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);

% the data.(odor) contains the binned spike counts (100 ms bins) for the
% different odors and trials (dimensions: unit x time x trials)

%
shuffle_labels = false;
[confusion_matrix, correct_or_not] = odor_response_classifier(pop_vec,shuffle_labels,zscore_data,ops);
%% compare to classifier trained with shuffled labels
% if compare_classifier
%     num_it = 3;
%    for i = 1:num_it
%        shuffle_labels = true;
%        [~, correct_or_not_shuffle] = odor_response_classifier(pop_vec,shuffle_labels,zscore_data,ops);
%        contingency_table = table([nnz(correct_or_not);nnz(correct_or_not_shuffle)],...
%            [nnz(~correct_or_not);nnz(~correct_or_not_shuffle)],'VariableNames',{'correct','not'},'RowNames',{'true','shuffle'});
%        [h,p(i),stats] = fishertest(contingency_table);
%    end
% end
shuffle_labels = true;
[~, correct_or_not_shuffle] = odor_response_classifier(pop_vec,shuffle_labels,zscore_data,ops);
contingency_table = table([nnz(correct_or_not);nnz(correct_or_not_shuffle)],...
   [nnz(~correct_or_not);nnz(~correct_or_not_shuffle)],'VariableNames',{'correct','not'},'RowNames',{'true','shuffle'});
[h,p,stats] = fishertest(contingency_table);
%% plot
f = figure;
if size(confusion_matrix,1)>3
    odor_names = ["ylang ylang"  "peanut butter" "BL6 #1" "BL6 #2" "CD1"];
else
    odor_names = ["fam" "nov1" "nov2"];
end
%%
curr_odors = odor_names(1:size(confusion_matrix,1));
heatmap(curr_odors,curr_odors,confusion_matrix,'GridVisible','off');
% imagesc(confusion_matrix);
ylabel('actual');
xlabel('predicted');
% title(['p:' num2str(median(p))]);
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

f.Units = 'centimeters';
% f.Position = [3 3 8.3 5.6];




