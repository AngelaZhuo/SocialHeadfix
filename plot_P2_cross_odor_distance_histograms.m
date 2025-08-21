function ff = plot_P2_cross_odor_distance_histograms(d,unit_maps,ops)
%%
% plot intra odor distance histogram as well
showIntra = 1; % plot intra odor distance histogram
zscored = 0;

%% Compute


% build spike count population vector
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);

% get cross-odor distance
[euclidian, cosine] = compute_popVec_distance_between_condition(pop_vec,ops);
if ops.basenorm_crossdist
    tmp = ops.response; ops.response = [0-diff(tmp) 0];
    [euclidian_base, cosine_base] = compute_popVec_distance_between_condition(pop_vec,ops);
    ops.response = tmp; clear tmp; 

    % element-wise division by cross trial distance of baseline
    euclidian = gdivide(euclidian,euclidian_base);
    cosine = gdivide(cosine,cosine_base);
end

% prepare timebase
% time_base = -ops.pre:ops.binsize:ops.post;
% time_base(end) = [];

% response and baseline bins
% response_bins = find(round(time_base,1)==ops.response(1)):find(round(time_base,1)==ops.response(2))-1;
% baseline_bins = find(round(time_base,1)==ops.baseline(1)):find(round(time_base,1)==ops.baseline(2))-1;

events = d.events{1,(d.map(unit_maps(1)))};
% odor_nums = unique([events.case_num]);
odor_num = unique([events([events.odor_dur]~=0).case_num]); %omit non-odor trials
if size(euclidian,1)<=3
    odor_names = ["f"  "n1" "n2"];
else
    odor_names = ["ylang ylang"  "peanut butter" "BL6 #1" "BL6 #2" "CD1"];
end
%% Plot
% plot estimate kernel density function for nicer visualization: ksdensity


f = figure;
f.Units = 'centimeters';

Y = nan(numel(euclidian{1,1}),sum((1:size(euclidian,1))-~showIntra));
counter = 1;
hold on
leg_strs= [];
for ii = 1:size(euclidian,1)
    for jj = 1:size(euclidian,2)
        if any([jj-showIntra<ii*zscored ii-showIntra<jj*~zscored])
%         if ii-showIntra < jj
            if ii ~= jj
                curr_eucl = reshape(euclidian{ii,jj},[],1);
                Y(:,counter) = curr_eucl;                
            else
                if zscored
                    curr_eucl = reshape(tril(euclidian{ii,jj})+tril(euclidian{ii,jj})',[],1);
                else
                    curr_eucl = reshape(triu(euclidian{ii,jj})+triu(euclidian{ii,jj})',[],1);
                end
                % curr_eucl = reshape(triu(euclidian{ii,jj}),[],1);
                curr_eucl(curr_eucl==0) = NaN;
                Y(:,counter) = curr_eucl;
            end
            counter = counter +1;
            if zscored || ops.basenorm_crossdist
                histogram(curr_eucl,'BinWidth',.1, 'Normalization', 'probability');
            else
                histogram(curr_eucl,'BinWidth',.5, 'Normalization', 'probability');
            end
            leg_strs = [leg_strs strjoin(sort([odor_names([ii  jj])])," vs ")];
        else
            continue;
        end
    end
end
leg=legend(leg_strs);
Y(:,all(isnan(Y),1)) = [];

%%
[p,tbl,stats] = anova1(Y,leg_strs);
figure
[c,m,h] = multcompare(stats);

fig_hds = findobj('type','figure');
if numel(fig_hds)>4
    fig_hds(5:end) = [];
end
for fx = 1:numel(fig_hds)
    fig_hds(fx).Position(1:2) = [3+10*(fx-1) 3];
end
ff = univCombFig(flip(fig_hds([1 2 4])),[1 3],ops.figvis,1);
ff.Children(1).Title.String='';
ff.Children(1).XLabel.String='';
close(fig_hds(:))
end