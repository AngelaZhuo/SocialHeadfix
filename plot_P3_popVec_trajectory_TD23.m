function f = plot_P3_popVec_trajectory_TD23(d,unit_maps, ops)
%% plots population vector trajectories over time and orthogonalizes for social vs. non-social odors
% build spike count population vector
% Modified from plot_P3_popVec_trajectory.m (DW)
[pop_vec,ops] = build_PopVec(d,unit_maps,ops);
% lineProps.col = {[1 0 1];[0.1 1 1];[1 1 0]};
% lineProps.col =     {[0 0 0.5156];[0.5000 1 0.5000];[0.5000 0 0]};
lineProps.width = 1;
plot3d= 1;
num_pcs = 3;
zscored = 1;

% trial_blocks = {1:29};
trial_blocks = {1:size(pop_vec{1,1},3)};
% trial_blocks = {1:9 10:18 19:27};
block_num = numel(trial_blocks);

if size(pop_vec,2)<=3
    odor_names = ["fam" "nov1" "nov2"];
    lineProps.col = get_colors(3);
    lineProps.col = cell2mat(repelem(lineProps.col,block_num)).*repmat((.5:.5/(block_num-1):1)',block_num,1);
else
    tmp=jet(256);
    lineProps.col = num2cell(tmp(1:floor(256/size(pop_vec,2)-1):256,:),2);%get_colors(1); % (non)social   
    lineProps.col = cell2mat(repelem(lineProps.col,block_num)).*repmat((.5:.5/(block_num-1):1)',block_num,1);
    odor_names = ["ylang ylang"  "peanut butter" "BL6 #1" "BL6 #2" "CD1"];
end
clear population averaged_population xe X mappedX
if zscored
    pop_vec(1,:) = [];
else
    pop_vec(2,:) = [];
end
population = cellfun(@(x) permute(x,[3 1 2]), pop_vec, 'UniformOutput', false);%cell for every conditions: trials x units x time
% average population across trials -> cell: units x time
for tb = 1:block_num
    for ph=1:size(population,2)
    % averaged_population{ph}(:,:)=mean(population{ph}(fvOn_bin_position:end,:,:),1);
    averaged_population{tb,ph}(:,:)=mean(population{ph}(trial_blocks{tb},:,:),1);
    end
end
if 1
    % temporal embedding
    m=4;
    tau=1;
    for ph=1:numel(averaged_population)
        if ~isempty(averaged_population{ph})
            xe{ph}=embedding(averaged_population{ph}',m,tau);
            xe{ph}=xe{ph}';
        end
    end
    X=[];
    for ph=1:numel(xe)
        X=[X,xe{ph}];
    end
else
    m = 0;
    X = [averaged_population{:}];
end
X(any(isnan(X),2),:)=[];
% dimensionality reduction: PCA
[coeff,score,latent] = pca(X');
mappedX=score(:,1:num_pcs);
%% Plot
f=figure;
% n_bins=size(xe{1},2);
n_bins = numel(-ops.pre:ops.binsize:ops.post)-max(m,1);
fvOn_bin_position = find(-ops.pre:ops.binsize:ops.post==0)-m;
fvOff_bin_position = find(-ops.pre:ops.binsize:ops.post==1)-m;
hold on;
start_bin = 1;
for ii = 1:size(population,2)*block_num
    if plot3d
        h(ii)=plot3(mappedX(start_bin:start_bin+n_bins-1,1), mappedX(start_bin:start_bin+n_bins-1,2), mappedX(start_bin:start_bin+n_bins-1,3), ...
            'Color', lineProps.col(ii,:), 'LineWidth', 1.5);
        
        % Plot start (point)
        plot3(mappedX(start_bin,1), mappedX(start_bin,2),mappedX(start_bin,3),'k.');
        % Plot fv_on (triangle)
        plot3(mappedX(start_bin+fvOn_bin_position-1,1), mappedX(start_bin+fvOn_bin_position-1,2),mappedX(start_bin+fvOn_bin_position-1,3),'k^');
        % Plot fv_off (square)
        plot3(mappedX(start_bin+fvOff_bin_position-1,1), mappedX(start_bin+fvOff_bin_position-1,2),mappedX(start_bin+fvOff_bin_position-1,3),'ks');
        % Plot end (circle)
        plot3(mappedX(start_bin+n_bins-1,1), mappedX(start_bin+n_bins-1,2),mappedX(start_bin+n_bins-1,3),'ko');
        
    else
        if ii == 1; tiledlayout(ceil(num_pcs/5),min([num_pcs 5])); end
        for pcx = 1:num_pcs
            nexttile(pcx); hold on
            h(ii) = plot(mappedX(start_bin:start_bin+n_bins-1,pcx),...
                'Color', lineProps.col(ii,:), 'LineWidth', 1.5);

            if ii==1
                try; xline(fvOn_bin_position); end
                xline(fvOff_bin_position)
                title(['PC' num2str(pcx)])
            elseif ii == 3                
            end
        end
    end
    start_bin = start_bin+n_bins;
end

hold off

if plot3d
    drawnow   
    grid on;
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
end
% set_fonts();
drawnow
legend(h(block_num:block_num:end),odor_names(1:size(population,2)),'FontSize',10,'LineWidth',2);
f.Units = 'centimeters';
% f.Position = [3 3 6 3];
end

%% subfunctions
function xe = embedding(averaged_population,m,tau)
xe = averaged_population(1+(m-1)*tau:end,:);
for mx = 1:m-1
    tmp = averaged_population(1+(m-1-mx)*tau:end-mx*tau,:);
    xe = cat(2,xe,tmp);
end
end