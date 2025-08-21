function f1 = plot_popVec_between_conditions(d,unit_maps,ops)
%%
odor_labels = ["fam" "nov1" "nov2"];
plot_norm =0;

%% Compute image


% build spike count population vector
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);

% calculate distance between conditions over bins
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
time_base = -ops.pre:ops.binsize:ops.post;
time_base(end) = [];

events = d.events{1,(d.map(unit_maps(1)))};
odor_nums = unique([events.case_num]);

%% Plot

f1 = figure('Position',[1 1 10 5*(1+plot_norm)]  ,'name','Population Vector between conditions');
% tiledlayout(1,2)
% fullscreen(f1);

% euclidian
s(1) = subplot(1+plot_norm,2,1);

temp = cell2mat(euclidian);
image_e = triu(temp)+triu(temp)';
imagesc(image_e);
colorbar;
title('Euclidian Distance','FontSize',10)
xy_lines_ticks(euclidian,odor_labels)

% cosine
s(2) = subplot(1+plot_norm,2,2+plot_norm);
temp = cell2mat(cosine);
image_c = triu(temp)+triu(temp)';
imagesc(image_c);
colorbar;
title('Cosine Distance','FontSize',10)
xy_lines_ticks(euclidian,odor_labels)

if plot_norm
    % euclidian norm
    s(3) = subplot(2,2,2);
    temp = cell2mat(euclidian);
    image_en = tril(temp)+tril(temp)';
    imagesc(image_en);
    colorbar;
    title('Normalized Euclidian Distance','FontSize',10)
    xy_lines_ticks(euclidian,odor_labels)

    % cosine normalized
    s(4) = subplot(2,2,4);
    temp = cell2mat(cosine);
    image_cn = triu(temp)+triu(temp)';
    imagesc(image_cn);
    colorbar;
    title('Normalized Cosine Distance','FontSize',10)
    xy_lines_ticks(euclidian,odor_labels)
end

% make pretty


% set(s,'XTick',10:20:90);set(s,'YTick',10:20:90);
% labels = cellfun(@get_event_label, num2cell(odor_nums), 'UniformOutput', false);
% 
% set(s,'XTicklabels',labels);
% set(s,'Yticklabels',labels);
set(s,'FontSize',8);

add_title(f1,[num2str(numel(unit_maps)),' units'],'FontSize',10);
% add_title(f1,['Population Vector between conditions (n=',num2str(numel(unit_maps)),')'],'FontSize',1);

end

%% subfunctions
function xy_lines_ticks(euclidian,labels)

lpos = cumsum(repmat(size(euclidian{1,1},1),1,size(euclidian,1)));
for lx = 1:numel(lpos)-1
    xline(lpos(lx)+.5)
    yline(lpos(lx)+.5)
end

tcks = lpos-diff(lpos(1:2))/2;
labels = labels(1:numel(tcks));
xticks(tcks); yticks(tcks);
xticklabels(labels); yticklabels(labels);
end