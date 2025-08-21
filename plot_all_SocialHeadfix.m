function f = plot_all_SocialHeadfix(d,uids,unit_labels,para_name,~,ops)
% modified by AZ 2025.08

svplots = 0;
svfgs = 0;
hmp = 0;
hst = 0;
bld = 0;
pct = 0;
cls = 1;

tx = ops.exclude_trials;
for i = 1:numel(uids)
    if ~isempty(unit_labels)
        t_string = {[para_name ' - ' unit_labels{i} ' - first ' num2str(numel(tx)) ' trials X - response ' num2str(ops.response(1)) ':' num2str(ops.response(2)) 'sec'];''};
        f_string = replace(replace(t_string{1},':','-'),' ','');
    end
    

    if hmp % cross odor heatmaps
        f = plot_popVec_between_conditions(d,uids{i},ops);
        ScaleFigure;drawnow;
        if ~isempty(unit_labels);add_title(f,t_string,'FontSize',12);end
        if svplots; exportgraphics(f,[save_dir 'euclidian_cosine_heatmaps_norm_' f_string '.png']); close(f); end
        if svfgs; savefig(f,[save_dir 'euclidian_cosine_heatmaps_norm_' f_string '.fig']); close(f); end
    end

    if hst % cross odor histograms
        plot_P2_cross_odor_distance_histograms(d,uids{i},ops);
        drawnow;f=gcf;
        if ~isempty(unit_labels);add_title(f,t_string,'FontSize',12);end
        if svplots; exportgraphics(f,[save_dir 'euclidian_hist_' f_string '.png']); close(f); end
        if svfgs; savefig(f,[save_dir 'euclidian_hist_' f_string '.fig']); close(f); end
    end

    if bld % euclidean from baseline
        f = plot_and_test_euclidean_from_base(d,uids{i},ops);
        if ~isempty(unit_labels);add_title(f,t_string,'FontSize',12);drawnow;end
        if svplots; exportgraphics(f,[save_dir 'euclidian_from_base_' f_string '.png']); close(f); end
        if svfgs; savefig(f,[save_dir 'euclidian_from_base_' f_string '.fig']); close(f); end
    end

    if pct % PC trajectory
        f = plot_P3_popVec_trajectory_TD23(d,uids{i}, ops); drawnow
        if ~isempty(unit_labels);add_title(f,[para_name ' - ' unit_labels{i}]);end
        if svplots; keyboard; exportgraphics(f,[save_dir 'PCtrajectory_' f_string '.png']); close(f); end
        if svfgs; savefig(f,[save_dir 'PCtrajectory_' f_string '.fig']); close(f); end
    end

    if cls % classifier
        [f,p] = plot_confusion_matrix_fam_nov(d,uids{i},ops);
        if ~isempty(unit_labels); title(['p:', num2str(p),' ', para_name, '-', unit_labels{i}]); end%add_title(f,[para_name ' - ' unit_labels{i}]);end
        if svplots; exportgraphics(f,[save_dir 'ConfusionMatrix_' f_string '.png']); close(f); end
        if svfgs; savefig(f,[save_dir  'ConfusionMatrix_' f_string '.fig']); close(f); end
    end
%     if ~svplots; keyboard; end

end
end