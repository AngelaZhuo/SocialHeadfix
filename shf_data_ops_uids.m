%% Script to load and flexibly switch between datasets
% modified by AZ 2025.08
%% set parameters
function [ops,d,para_name,unit_labels,uids]= shf_data_ops_uids(setops,d_select)
if ~exist('setops','var');if ~exist('ops','var');setops=1;else;setops=input('Reset ops? 0 or 1 \n');end; end
if setops
    ops.pre = 1;
    ops.post = 4;
    ops.binsize = .05;
    ops.baseline = [-1 0];
    ops.response =[0.05 1];
    ops.basenorm_crossdist = 1;
    ops.sniff_align = 0;
    ops.average_trials = 1;
    ops.smoothing = 0;
    ops.shuffle_trial_order=0;
    ops.align_to_event='fv_on';
    ops.event_name ='fv_on';
    % ops.exclude_trials =1:3;
    ops.exclude_trials = [];
    ops.shuffle_conditions = 0;
    ops.figvis = 1;
end
clear setops

%% load or select dataset(s) and the the respective unit vectors
if ~exist('d_select','var')
    d_select = {input([...
    'To load all datasets press enter.\n\n' 'To get a specific dataset enter one of the following strings, then press enter:\n' ...
    'OXT_1fam1nov\n' ...
    'OXT_(non)social\n' ...
    'TD23_1fam2nov\n' ...
    'OFC-AI_1fam1nov\n' ...
    'OFC-AI_1fam2nov\n' ...
    'PPC-LEC_1fam1nov\n' ...
    'PPC-LEC_(non)social\n' ...
    'reappraisal\n'...
    'social_defeat\n'...
    ])};
end
if ~exist('d','var'); d.info.name = ''; end

if isempty(d_select{1})
    d_select = {'OXT_1fam1nov' ...
        'OXT_(non)social' ...
        'TD23_1fam2nov' ...
        'OFC-AI_1fam1nov' ...
        'OFC-AI_1fam2nov' ...
        'PPC-LEC_1fam1nov' ...
        'PPC-LEC_(non)social'...
        'reappraisal'...
        'social_defeat'...
        };
end

for dx = 1:numel(d_select)
    current_dataset = d_select{dx};
    switch current_dataset
        case 'OXT_1fam1nov'
            %% OXT

            % % groups
            % % •	Gruppe 1 (AON optogenetic OXT boost): group==1 & ref_viols<2 & mean_fr > 0.1 & mean_fr < 20 & ~exclude
            % % •	Gruppe 2 (AON control): group==2 & ref_viols<2 & mean_fr > 0.1 & mean_fr < 20 & ~exclude
            % % •	Gruppe 3 (4 in unit_maps) (AON OXTR knockout): group==3 & ref_viols<2 & mean_fr > 0.1 & mean_fr < 20 & ~exclude
            % % •	Gruppe 4 (3 in unit_maps) (AON OXTR knockout control == dTom injection anstatt Cre): group==4 & ref_viols<2 & mean_fr > 0.1 & mean_fr < 20 & ~exclude
            % % •	VTA: trode >= 33 & ref_viols<2 & mean_fr > 0.1 & mean_fr <= 12 & ~exclude
            % % •	VS: trode < 33 & ref_viols<2 & mean_fr > 0.1 &  ~exclude
            % %
            if ~exist('d_1f1n_oxt','var')
                if ~all(ismember(unique({d.info.name}),{'OXT10','OXT16','OXT23','OXT28','oxt10','oxt12','oxt16','oxt17', ...
                        'oxt21','oxt23','oxt24','oxt25','oxt26','oxt27','oxt28','oxt30','oxt31','oxt32','oxt33','oxt34', ...
                        'srm010','srm011','srm012','srm013','srm014','srm015','srm016','srm017','srm018','srm019','srm020', ...
                        'srm021','tdat01','tdat02','tdat04'}))
                    disp('Loading data:')
                    load('/zi/flstorage/group_entwbio/data/Mirko/SocialHeadFix/1Fam1Novel/OXT/250312_data.mat');
                end
                d_1f1n_oxt = d;
            else
                if ~all(ismember(unique({d.info.name}),{'OXT10','OXT16','OXT23','OXT28','oxt10','oxt12','oxt16','oxt17', ...
                        'oxt21','oxt23','oxt24','oxt25','oxt26','oxt27','oxt28','oxt30','oxt31','oxt32','oxt33','oxt34', ...
                        'srm010','srm011','srm012','srm013','srm014','srm015','srm016','srm017','srm018','srm019','srm020', ...
                        'srm021','tdat01','tdat02','tdat04'}))
                    d = d_1f1n_oxt;
                end
            end
            unit_maps =load('/zi/flstorage/group_entwbio/data/Mirko/SocialHeadFix/1Fam1Novel/OXT/unit_maps_OXT2024.mat'); % ref_viol and frequency corrected;
            unit_labels_oxt = fieldnames(unit_maps);
            unit_labels = cell(1);
            uids = cell(1);
            % for ulx = 1:numel(unit_labels_oxt)
            for ulx = 1:numel(unit_labels_oxt)-2
                uids{ulx} = unit_maps.(unit_labels_oxt{ulx});
                unit_labels{ulx} = replace(unit_labels_oxt{ulx},'unit_maps_','');
            end
            unit_labels{ulx+1} = 'pMSN';
            uids{ulx+1} = find(contains({d.clust_params.animal},'tdat') & ...
                [d.clust_params.mean_fr]>.1 & [d.clust_params.mean_fr]<=5 & ...
                [d.clust_params.trode]<=32 & [d.clust_params.ref_viols]<2);

            unit_labels{ulx+2} = 'pFSI';
            uids{ulx+2} = find(contains({d.clust_params.animal},'tdat') & ...
                [d.clust_params.mean_fr]>5 & [d.clust_params.trode]<=32 & ...
                [d.clust_params.ref_viols]<2);

            unit_labels{ulx+3} = 'VTA(pDAN)';
            uids{ulx+3} = find(contains({d.clust_params.animal},'tdat') & ...
                [d.clust_params.mean_fr]>1 & [d.clust_params.mean_fr]<=12 & ...
                [d.clust_params.trode]>32 & [d.clust_params.ref_viols]<2& ...
                ~ismember(1:numel(d.spikes),3844));

            unit_labels{ulx+4} = 'VTA(fast)';
            uids{ulx+4} = find(contains({d.clust_params.animal},'tdat') & ...
                [d.clust_params.mean_fr]>12 & [d.clust_params.trode]>32 & ...
                [d.clust_params.ref_viols]<2);
            uids = uids';
            para_name = '1fam1nov'
            unit_labels=unit_labels'
            % save_dir = '1Fam1Novel/OXT/'
            % if ~isfolder(save_dir); mkdir(save_dir);end
            clear unit_maps unit_labels_oxt

        case 'OXT_(non)social'
            %% OXT social nonsocial
            % load '/zi/flstorage/group_entwbio/data/David/data/OXT DATA 2019/P2 social vs non-social/data/data.mat';

            if ~exist('d_s_ns_oxt','var')
                if ~all(ismember(unique({d.info.name}),{'oxt10' 'oxt12' 'oxt16' 'oxt17' 'oxt21' 'oxt23' 'oxt24' 'oxt25' 'oxt26' 'oxt27' 'oxt28' 'oxt30'...
                        'oxt31' 'oxt32' 'oxt33' 'oxt34' 'srm010' 'srm011' 'srm012' 'srm013' 'srm014' 'srm015' 'srm016' 'srm017'...
                        'srm018' 'srm019' 'srm020' 'srm021' 'tdat01' 'tdat02' 'tdat04' 'tdat05' 'tdat06'}))
                    disp('Loading data:')
                    load('/zi/flstorage/group_entwbio/data/David/data/OXT DATA 2019/P2 social vs non-social/data/data.mat');
                end
                d_s_ns_oxt = d;
            else
                if ~all(ismember(unique({d.info.name}),{'oxt10' 'oxt12' 'oxt16' 'oxt17' 'oxt21' 'oxt23' 'oxt24' 'oxt25' 'oxt26' 'oxt27' 'oxt28' 'oxt30'...
                        'oxt31' 'oxt32' 'oxt33' 'oxt34' 'srm010' 'srm011' 'srm012' 'srm013' 'srm014' 'srm015' 'srm016' 'srm017'...
                        'srm018' 'srm019' 'srm020' 'srm021' 'tdat01' 'tdat02' 'tdat04' 'tdat05' 'tdat06'}))
                    d = d_s_ns_oxt;
                end
            end
            clear uids unit_labels

            i = 0;
            i=i+1; unit_labels{i} = 'AON';
            uids{i} =find(contains({d.clust_params.animal},'oxt','IgnoreCase',true) &...
                [d.clust_params.mean_fr]>.1 & [d.clust_params.mean_fr]<=20 & ...
                [d.clust_params.trode]<17 & [d.clust_params.ref_viols]<2 & ...
                ~ismember(d.map,[13])); %% no fv_on time in events
            i=i+1; unit_labels{i} = 'pDAN';
            uids{i} =find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>1 & [d.clust_params.mean_fr]<=12 & ...
                [d.clust_params.trode]>32 & [d.clust_params.ref_viols]<2);
            i=i+1; unit_labels{i}=  'VTAfast';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>12 &...
                [d.clust_params.trode]>32 & [d.clust_params.ref_viols]<2);
            i=i+1; unit_labels{i} = 'NAcc(pMSN)';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>.1 & [d.clust_params.mean_fr]<5  & ...
                [d.clust_params.ref_viols]<2 & [d.clust_params.trode]<17);
            i=i+1; unit_labels{i} = 'Tu(pMSN)';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>.1 & [d.clust_params.mean_fr]<5  & ...
                [d.clust_params.ref_viols]<2 &...
                [d.clust_params.trode]>16 & [d.clust_params.trode]<33);
            i=i+1; unit_labels{i} = 'NAcc(fast)';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>5  & [d.clust_params.ref_viols]<2 &...
                [d.clust_params.trode]<17);
            i=i+1; unit_labels{i} = 'Tu(fast)';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>5  & [d.clust_params.ref_viols]<2 &...
                [d.clust_params.trode]>16 & [d.clust_params.trode]<33);
            i=i+1; unit_labels{i} = 'VS(fast)';
            uids{i}=find(contains({d.clust_params.animal},'tdat') &...
                [d.clust_params.mean_fr]>5  & [d.clust_params.ref_viols]<2 &...
                [d.clust_params.trode]<33);
            i=i+1; unit_labels{i} = 'pDAN(an124)';
            uids{i} = uids{1}(ismember({d.clust_params(uids{strcmp(unit_labels,'pDAN')}).animal},{'tdat01', 'tdat02','tdat04'}));
            i=i+1; unit_labels{i} = 'VTAfast(an124)';
            uids{i} = uids{2}(ismember({d.clust_params(uids{strcmp(unit_labels,'VTAfast')}).animal},{'tdat01', 'tdat02','tdat04'}));
            uids = reshape(uids,[],1);

            para_name = '(non)social'
            unit_labels= reshape(unit_labels,[],1)
            % save_dir = '(non)social/OXT/'
            % if ~isfolder(save_dir); mkdir(save_dir);end
            

        case 'TD23_1fam2nov'
            %% TD23 (1fam2novel)
            if ~exist('d_1f2n_vs_vta','var')
                if ~all(ismember(unique({d.info.name}),{'x01','x02','x03','x04','x05','x06','x07','x08','x09','x10'}))
                    disp('Loading data:')
                    load('/zi/flstorage/group_entwbio/data/Angela/DATA/TD23/D-struct/KS3/social/d_social_11-Mar-2025.mat');
                end
                d_1f2n_vs_vta = d;
            else
                d = d_1f2n_vs_vta;
            end

            para_name = '1fam2nov'
            unit_labels= {'NAcc(pMSN)' 'NAcc(pFSI)' 'Tu(pMSN)' 'Tu(pFSI)' 'VS-VP fast' 'pDAN' 'VTAfast'}'
            % save_dir = '1Fam2Novel/VS-VTA/'
            % if ~isfolder(save_dir); mkdir(save_dir);end
            uids = cell(1);
            uids{1} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>.25 & ...
                [d.clust_params.mean_fr]<=5 & [d.clust_params.region_coding]==1);
            uids{2} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>5 & ...
                [d.clust_params.region_coding]==1);
            uids{3} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>.25 & ...
                [d.clust_params.mean_fr]<=5 & [d.clust_params.region_coding]==2 & ...
                ~ismember(1:numel(d.spikes),[695 741])); %exclude outlier units 
            uids{4} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>5 & ...
                [d.clust_params.region_coding]==2);
            uids{5} = [uids{2} uids{4}];  
            uids{6} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>1 & ...
                [d.clust_params.mean_fr]<=12 & [d.clust_params.region_coding]==3 & ...
                ~ismember(1:numel(d.spikes),[102 530]));
            uids{7} = find([d.clust_params.ref_viols]<2 & [d.clust_params.mean_fr]>12 & ...
                [d.clust_params.region_coding]==3 &  ~ismember(1:numel(d.spikes),[198 585 768]));
            uids = uids';

        case 'OFC-AI_1fam1nov'
            %% OFC AI 1f1n
            if ~exist('d_1f1n_ofc_ai','var')
                if ~all(ismember(unique({d.info.name}), {'d01','d02','d04','d05','d06','d07','d08','d09','d10'}))||...
                        numel(unique([d.events{1}([d.events{1}.odor_dur]~=0).case_num]))~=2
                    disp('Loading data:')
                    load('/zi/flstorage/group_entwbio/data/Mirko/SocialHeadFix/1Fam1Novel/OFC-AI/DATA/data_250409.mat')
                end
                d_1f1n_ofc_ai = d;
            else
                if ~all(ismember(unique({d.info.name}), {'d01','d02','d04','d05','d06','d07','d08','d09','d10'}))||...
                        numel(unique([d.events{1}([d.events{1}.odor_dur]~=0).case_num]))~=2
                    d = d_1f1n_ofc_ai;
                end
            end

            para_name = '1fam1nov'
            unit_labels = {'OFC' 'AI'}'
            % save_dir = '1Fam1Novel/OFC-AI/'
            % if ~isfolder(save_dir); mkdir(save_dir);end
            uids = cell(1);
            uids{1} = find([d.clust_params.mean_fr]>.25 & ...
                ismember([d.clust_params.trode], [8,9,10,12,13,15,16]) &...
                [d.clust_params.ref_viols]<2);
            uids{2} = find([d.clust_params.mean_fr]>.25 & ...
                ismember([d.clust_params.trode],  [1,2,3,4,5,6,7,11,14]) &...
                [d.clust_params.ref_viols]<2);
            
        case 'OFC-AI_1fam2nov'
            %% OFC AI 1f2n
            if ~exist('d_1f2n_ofc_ai','var')
                if ~all(ismember(unique({d.info.name}), {'d01','d02','d04','d05','d06','d07','d08','d09','d10'}))||...
                        numel(unique([d.events{1}([d.events{1}.odor_dur]~=0).case_num]))~=3
                    disp('Loading data:')
                    load('/zi/flstorage/group_entwbio/data/Mirko/SocialHeadFix/1Fam2Novel/OFC-AI/DATA/data_250409.mat')
                end
                d_1f2n_ofc_ai = d;
            else
                if ~all(ismember(unique({d.info.name}), {'d01','d02','d04','d05','d06','d07','d08','d09','d10'}))||...
                        numel(unique([d.events{1}([d.events{1}.odor_dur]~=0).case_num]))~=3
                    d = d_1f2n_ofc_ai;
                end
            end

            para_name = '1fam2nov'
            unit_labels = {'OFC' 'AI'}'
            % save_dir = '1Fam2Novel/OFC-AI/'
            % if ~isfolder(save_dir); mkdir(save_dir);end
            uids = cell(1);
            uids{1} = find([d.clust_params.mean_fr]>.25 & ...
                ismember([d.clust_params.trode], [8,9,10,12,13,15,16]) &...
                [d.clust_params.ref_viols]<2);
            uids{2} = find([d.clust_params.mean_fr]>.25 & ...
                ismember([d.clust_params.trode],  [1,2,3,4,5,6,7,11,14]) &...
                [d.clust_params.ref_viols]<2);
            
        case 'PPC-LEC_1fam1nov'
            %% ppc lec
            if ~exist('d_1f1n_ppc_lec','var')
                if    ~all(ismember(unique({d.info.name}),{'Z16','Z07','Z08','Z10','Z11','Z11','Z13','Z14','Z15','Z17'}))
                    disp('Loading data:')
                    load('//zi/flstorage/group_entwbio/data/Angela/DATA/z-cohort (PPC&LEC)/PROCESSED/Processed Dataperirhinal/P3/data_250311.mat','d')
                end
                d_1f1n_ppc_lec = d;
            else
                if    ~all(ismember(unique({d.info.name}),{'Z16','Z07','Z08','Z10','Z11','Z11','Z13','Z14','Z15','Z17'}))
                    d = d_1f1n_ppc_lec;
                end
            end

            para_name = '1fam1nov'
            unit_labels = {'PPC' 'LEC'}'
            % save_dir = '1Fam1Novel/PPC-LEC/'
            % if ~isfolder(save_dir); mkdir(save_dir);end

            uids = cell(1);
            uids{1} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.exclude_unit] == 0 & ...
                [d.clust_params.region] == 1);

            uids{2} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.exclude_unit] == 0 & ...
                [d.clust_params.region] == 2);
            

        case 'PPC-LEC_(non)social'
             %% ppc lec (non)social
            if ~exist('d_s_ns_ppc_lec','var')
                if    ~all(ismember(unique({d.info.name}),{'Z08' 'Z09' 'Z10' 'Z11' 'Z12' 'Z13' 'Z14' 'Z15' 'Z16' 'Z17'}))
                    disp('Loading data:')
                    load('//zi/flstorage/group_entwbio/data/Angela/DATA/z-cohort (PPC&LEC)/PROCESSED/Processed Dataperirhinal/P2/data_250408.mat')
                end
                d_s_ns_ppc_lec = d;
            else
                if    ~all(ismember(unique({d.info.name}),{'Z08' 'Z09' 'Z10' 'Z11' 'Z12' 'Z13' 'Z14' 'Z15' 'Z16' 'Z17'}))
                    d = d_s_ns_ppc_lec;
                end
            end

            para_name = '(non)social'
            unit_labels = {'PPC' 'LEC'}'
            % save_dir = '(non)social/PPC-LEC/'
            % if ~isfolder(save_dir); mkdir(save_dir);end

            uids = cell(1);
            uids{1} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.region] == 1);

            uids{2} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.region] == 2);
       
        case 'social_defeat'
             %% OFC AI (danai) VS-VTA (Angela)
            if ~exist('d_sdf','var')
                if ~all(ismember(({d.info.paradigm}),'social_defeat'))
                    disp('Loading data:')
                    disp('file not sved yet')
                    % load('//zi/flstorage/group_entwbio/data/Angela/DATA/z-cohort (PPC&LEC)/PROCESSED/Processed Dataperirhinal/P2/data_250408.mat')
                end
                d_sdf = d;
            else
                if ~all(ismember(({d.info.paradigm}),'social_defeat'))
                    d = d_sdf;
                end
            end

            para_name = 'social_defeat'
            unit_labels = {'OFC' 'AI' 'NAcc(pMSN)' 'NAcc(pFSI)' 'Tu(pMSN)' 'Tu(pFSI)' 'VS-VP fast' 'pDAN' 'VTAfast'}
            % save_dir = '(non)social/PPC-LEC/'
            % if ~isfolder(save_dir); mkdir(save_dir);end

            uids = cell(1);
            uids{1} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.region] == 1);

            uids{2} = find([d.clust_params.ref_viols] < 2 & ...
                [d.clust_params.mean_fr] > .1 & ...
                [d.clust_params.mean_fr] < 30 & ...
                [d.clust_params.region] == 2);

        otherwise
            warning([current_dataset ' not implemented'])
    end

end
    clear d_select i dx ulx
end