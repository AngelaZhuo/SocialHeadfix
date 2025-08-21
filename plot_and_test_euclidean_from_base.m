function f = plot_and_test_euclidean_from_base(d,unit_maps,ops)
base_trialwise = 0;
%% compute pop_vec
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);
%% compute euclidian
euclidian = compute_popVec_distance_from_baseline(pop_vec,ops);

%% plot
f = figure('Position',[1 1 20 10]);
odor_names = ["fam" "nov1" "nov2"];
%% line plots
subplot(1,2,1)
% time_base = round(-ops.pre:ops.binsize:ops.post,-log10(ops.binsize));
time_base = round(-ops.pre:ops.binsize:ops.post,3);
e_idx = 1 + base_trialwise;
% l_colors = [1 0 0; 0 1 0; 0 0 1];
l_colors = cat(1,cell2mat(get_colors(3)));
for ii = 1:size(euclidian,2)
    M = mean(mean(euclidian{e_idx,ii},3),1);
    SE = std(mean(euclidian{e_idx,ii},3),0,1)/sqrt(size(euclidian{e_idx,ii},1));
    boundedline(time_base(1:end-1),M,SE*1.96,'cmap',l_colors(ii,:),'alpha')
end
ax = gca;
patch([ops.response flip(ops.response)],reshape([ylim;ylim],1,[]),[.1 .1 .1],'FaceAlpha',.25,'EdgeColor','none')
xline(0,'--');xline(1,'--')
legend(flip(ax.Children([3 4:2:end])),[odor_names(1:size(euclidian,2)) "test window"])

if base_trialwise
    title('Euclidian distance from baseline (trialwise)')
else
    title('Euclidian distance from baseline (average)')
end
%% boxplots
subplot(1,2,2)
% time_base = round(-ops.pre:ops.binsize:ops.post,-log10(ops.binsize));
time_base = round(-ops.pre:ops.binsize:ops.post,3);
response_bins = find(time_base==ops.response(1)):find(time_base==ops.response(2));


for ii = 1:size(euclidian,2)
    mean_euclidian{1,ii} = nanmean(euclidian{e_idx,ii}(:,response_bins),2);
end
p = NaN(1,sum((1:ii)-1));
i = 1;
pairs = [];
for ii = 1:size(mean_euclidian,2)
    for ij = 1:size(mean_euclidian,2)
        if ii < ij
            pairs = [pairs;ii ij];
            [~,p(i)]=ttest(mean_euclidian{ii},mean_euclidian{ij});
            i = i+1;
        end
    end
end

boxplot([mean_euclidian{:}])
xticklabels(odor_names)


yt = get(gca, 'YTick');
axis([xlim   ylim.*[1 1.3]])
xt = get(gca, 'XTick');
hold on
for px = 1:size(pairs,1)
    if p(px)<=.001
        marker_color = [0 0 0];
    elseif p(px)<=.05
        marker_color = [.5 .5 .5];
    else
        continue
    end
    plot(xt(pairs(px,:)), [1 1]*max(yt)*(1+.075*px), '-k', ...
        mean(xt(pairs(px,:))), max(yt)*(1.025+.075*px), '*', 'Color', marker_color)

end
hold off
title('grey: p<0.05     |     black: p<0.001')

end