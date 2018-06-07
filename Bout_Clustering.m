% Bout_Clustering 

% marcus.ghosh.11@ucl.ac.uk 

%% Info 

% Clusters both active and inactive bouts from viewpoint data, 
% then makes plots & calculates statistics relating to the clusters 

% Input 
    % Workspace or workspaces output from Vp_Extract.m  

% Output 
    % Figures & Statistics 
    % Workspace that can be used with x.m to look @ transition patterns
        % between clusters 

%% Assumptions 

%% Dependencies 

% Files  
    % Lbmap
        %https://uk.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps

    % Knee_pt
        %https://uk.mathworks.com/matlabcentral/fileexchange/35094-knee-point
        
% MATLAB Products 
    % Statistics and Machine Learning Toolbox
    
%% Options     
set(0,'DefaultFigureWindowStyle','docked'); % dock figures
set(0,'defaultfigurecolor',[1 1 1]); % white background

%% Selecting Files  

% Load Data - using multiselect
[filename, pathname] = uigetfile('*.mat', 'Select files','MultiSelect','on'); %Select files
if isequal(filename,0) %If no file is selected
    error('No Files Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
end

%% Load bout structure data (Parameter Extracted Data)

tic
% Data Structures (Concatenate Variables)
% Note for (2,1) cells - 1 = wake, 2 = sleep
wake_cells = [];
sleep_cells = [];
i_group_tags = [];
i_experiment_tags = [];
experiment_tags = cell(2,1);
group_tags = cell(2,1);
fish_tags = cell(2,1);
parameter_indicies = cell(2,1);

counter = 1; % Start a counter
for f = 1:size(filename,2) %For each file
    clear experiment;
    experiment = load(strcat(pathname,filename{f})); %Load the mat file
    
    % Nab variables
    days_crop{f} = experiment.days_crop; % days crop
    nights_crop{f} = experiment.nights_crop; % nights crop
    parameters{f} = experiment.parameters; % parameters
    cmap{f} = experiment.cmap; % color map
    cmap_2{f} = experiment.cmap_2; % expanded color map
    night_color{f} = experiment.night_color; % night color
    geno_list{f} = experiment.geno_list; % group names
    units{f} = experiment.units; % units
    unit_conversion{f} = experiment.unit_conversion; % unit conversion
    days{f} = experiment.days; % days
    nights{f} = experiment.nights; % nights
    first_night{f} = experiment.first_night; % first night
    time_window{f} = experiment.time_window; % time window 
    fps{f} = experiment.fps; % fps 
    lb{f} = experiment.lb; % lb
    lb_sec{f} = experiment.lb_sec; % lb in seconds 
    
    % Concatenate variables
    for i = 1:size(experiment.wake_cells,2) % For each fish
        wake_cells = [wake_cells ; experiment.wake_cells{1,i}]; % wake cells
        sleep_cells = [sleep_cells ; experiment.sleep_cells{1,i}]; % sleep cells
        parameter_indicies{1,1} = [parameter_indicies{1,1} ; experiment.parameter_indicies{1,i}]; % wake bout windows
        parameter_indicies{2,1} = [parameter_indicies{2,1} ; experiment.parameter_indicies{2,i}]; % sleep bout windows
        group_tags{1,1} = [group_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*...
            experiment.group_tags(i,1)]; % wake group tags
        group_tags{2,1} = [group_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*...
            experiment.group_tags(i,1)]; % sleep group tags
        experiment_tags{1,1} = [experiment_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*f]; % wake experiment tags
        experiment_tags{2,1} = [experiment_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*f]; % sleep experiment tags
        fish_tags{1,1} = [fish_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*counter]; % wake fish tags
        fish_tags{2,1} = [fish_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*counter]; % sleep fish tags
        counter = counter + 1; % add to fish counter
    end
    
    %delta_px_sq{1,f} = experiment.delta_px_sq; skipping this is easier on memory
    i_group_tags = [i_group_tags ; experiment.group_tags]; % individual fish group tags
    i_experiment_tags = [i_experiment_tags ; ones(size(experiment.wake_cells,2),1)*f]; % individual fish experiment tags
    
end

clear experiment f i counter;
toc

%% Options 

% Grouping repeats of experiments (hard coded)
experiment_reps = [1 1 1 2 2 3 4]; % experiment groupings 
i_experiment_reps = i_experiment_tags;
for er = 1:max(experiment_reps) % for each repeat 
    found = find(experiment_reps == er); % find experiments
    
    for f = found % for each experiment in the repeat 
        i_experiment_reps(i_experiment_reps == f,1) = er; % tag with grouping variable 
    end 
    
end 

% Adjust colours 
for e = 1:size(cmap,2) % for each experiment
    if max(i_group_tags(i_experiment_tags==e)) == 1 % if theres just one group (e.g. WT experiments)
        cmap{e}(1,:) = [135 206 250]/255; % light sky blue
        cmap_2{e}(1,:) = cmap{e};
        cmap_2{e}(2,:) = [25 25 112]/255; % midnight blue
    else % for experiments with multiple groups 
        cmap_2{e} = flip(cmap_2{e}); % flip cmap (so it starts with blue)
        for c = 1:2:size(cmap_2{e},1) % for every other color 
            cmap_2{e}([c c+1],:) = cmap_2{e}([c+1 c],:); % swap pairs of colours around 
        end
        cmap{e} = cmap_2{e}(1:2:size(cmap_2{e},1),:); % Extract main colors
    end
end

% Adjust time windows (Hard coded) 
for e = 1:size(time_window,2) % for each experiment 
    if experiment_reps(e) < 3 
       time_window{e} = [3 6]; % take the middle two days/nights 
       days{e} = [2 3]; nights{e} = [2 3]; 
    else
       time_window{e} = [1 2]; % take the first day/night 
       days{e} = 1;
    end 
end 

clear er found f e c 

%% Preparing Data for Clustering  

% Active 
tic
X{1,1} = []; % empty X
% Z-score each fishes data
for f = 1:max(fish_tags{1,1}) % for each fish
    X{1,1} = [X{1,1} ; zscore(wake_cells(fish_tags{1,1} == f,3:end))];
    if mod(f,100) == 0 % report every 100 fish 
        disp(horzcat('Completed ',num2str(f),' fish of ',...
            num2str(max(fish_tags{1,1}))));
    end
end
toc

[coeff,score,~,~,explained,~] = pca(X{1,1}); % pca
% Note: By default Matlab Centers data for PCA by subtracting
% the mean from each column (as the means are not quite zero this seems
% appropriate)
[knee_dim] = knee_pt(explained); % Choose this many dimensions
disp(horzcat('Reduced Wake data to ',num2str(knee_dim),' dimensions',...
    ' Explains ',num2str(sum(explained(1:knee_dim))),' % '));
X{1,1} = score(:,1:knee_dim);

% Inactive 
% Handling NaN Values 
    % Note that there are so few NaN values that giving them "temp" values 
    % for the clustering won't make a difference 
sleep_cells_nan_track = isnan(sleep_cells(:,3)); % store nan locations  
sleep_cells(sleep_cells_nan_track,3) = 1; % set NaN's to 1 (the mode of the data) 

X{2,1} = []; % empty X
X{2,1} = sleep_cells(:,3); 

clear f

%% PCA Figure 

figure; 
% Scree Plot 
subplot(1,2,1); hold on; 
plot(1:length(explained),explained,'linewidth',9,'color',night_color{1});
scatter(1:length(explained),explained,270,'filled',...
    'markerfacecolor',night_color{1},'markeredgecolor','k');
scatter(1,explained(1),270,'filled','markerfacecolor',cmap_2{1}(1,:)); 
scatter(2,explained(2),270,'filled','markerfacecolor',cmap_2{1}(2,:)); 
scatter(3,explained(3),270,'filled','markerfacecolor',[1 0.5 0]); 
xlabel('Component','Fontsize',32); 
ylabel('Variance Explained (%)','Fontsize',32);
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
set(gca,'FontName','Calibri');
axis([1 length(explained) ylim]);

% Coeff Plot 
subplot(1,2,2); hold on; 
plot([1,length(explained)],[0 0],'--k','linewidth',3); 
plot(coeff(:,1),'color',cmap_2{1}(1,:),'linewidth',3);
plot(coeff(:,2),'color',cmap_2{1}(2,:),'linewidth',3); 
plot(coeff(:,3),'color',[1 0.5 0],'linewidth',3); 
ylabel('Coefficient','Fontsize',32); 
set(gca,'XTick',1:length(explained)); 
set(gca,'XTickLabels',{'Length','Mean','Std','Total','Min','Max'},'fontsize',32); 
xtickangle(45); 
axis tight; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
set(gca,'FontName','Calibri');

% Save here 

%% Clustering Settings 
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat', 'X'); % load X
reps = 200; % set the number of repetitions
k_vals = 2:20; % set values of k (clusters) to try
a_size = 40000; % number of probe points  
s_vals = [40000,100000]; % min & max points to sample (uniformly)
GMM_reps = 5; % number of GMM Models to fit per iteration 
max_its = 1000; % number of GMM iterations (per repetition) 
method = 'average'; % linkage measure 
nn = 50; % number of nearest neighbours

% Calculate Regularization
score_values = unique(X{1,1}(:)')'; % Find unique values
score_zero = knnsearch(score_values,0); % Find the closest to zero
rv = abs(score_values(score_zero)); % Regularization value for GMM 

clear score_values score_zero 

%% Clustering 
save_pathname = uigetdir([],'Select a save location'); 

for s = 1:2 % for active & inactive
    tic 
    [ea, idx, idx_cts, ~, ...
        ea_links, ea_idx, ~, th, sample_a,sample_a_n] = ...
        gmm_sample_ea(X{s,1},reps,k_vals,a_size,s_vals,rv,GMM_reps,max_its,method,nn);
    toc 
    save(strcat(save_pathname,'\','180516','_',num2str(s),'.mat'),'-v7.3'); % save data  
    clear ea idx idx_cts ea_links ea_idx th sample_a sample_a_n
end

clear s 

%% Load Clustered Data  
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat'); 
active = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_1.mat',...
    'ea','idx','idx_cts','ea_links','ea_idx','th','sample_a','sample_a_n'); 
inactive = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_2.mat',...
    'ea','idx','idx_cts','ea_links','ea_idx','th','sample_a','sample_a_n');

% merge variables 
ea{1,1} = active.ea; ea{2,1} = inactive.ea; 
idx{1,1} = active.idx; idx{2,1} = inactive.idx; 
idx_cts{1,1} = active.idx_cts; idx_cts{2,1} = inactive.idx_cts; 
ea_links{1,1} = active.ea_links; ea_links{2,1} = inactive.ea_links; 
ea_idx{1,1} = active.ea_idx; ea_idx{2,1} = inactive.ea_idx; 
th(1,1) = active.th; th(2,1) = inactive.th; 
sample_a{1,1} = active.sample_a; sample_a{2,1} = inactive.sample_a; 
sample_a_n{1,1} = active.sample_a_n; sample_a_n{2,1} = inactive.sample_a_n; 

clear active inactive

%% Post-Clustering Settings 

% Colormap
numComp = [max(ea_idx{1,1}) max(ea_idx{2,1})]; % number of active & inactive numbers of clusters
scrap = lbmap(sum(numComp),'RedBlue'); % Color Scheme
cmap_cluster{1,1} = flip(scrap(1:numComp(1),:)); % generate a colormap for the clusters 
cmap_cluster{2,1} = scrap((numComp(1)+1):end,:); % generate a colormap for the clusters 
strings{1,1} = 'Active'; strings{2,1} = 'Inactive'; 

% Remove NaN Values from Sleep matracies  
sleep_cells(sleep_cells_nan_track,3) = NaN; % sleep cells 
idx{2,1}(sleep_cells_nan_track,1) = NaN; % cluster assignments  

% Merge Variables for ease of Looping 
cells{1,1} = wake_cells; 
cells{2,1} = sleep_cells; % cells 

clear scrap 

%% Sorting Clusters by Mean Length 
             
for s = 1:2 % for active & inactive
    mean_cluster_length = nan(1,numComp(s),'single'); % pre-allocate
    for c = 1:numComp(s) % For each cluster
        mean_cluster_length(c) = nanmean(cells{s,1}(idx{s,1}==c,3));
        % Calculate mean bout length
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    idx_numComp_sorted{s,1} = nan(size(cells{s,1},1),1,'single'); % Pre-allocate
    
    for c = 1:numComp(s)  % For each cluster
        idx_numComp_sorted{s,1}(idx{s,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    clear s mean_cluster_length c O;
end

clear idx

%% Evidence Accumulation Figure 
    % NOTE: Either need to edit colours to match numComp_sorted or go for a
    % different color scheme 
    
for s = 1:2 % for active & inactive
    figure; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    set(gca,'FontName','Calibri');
    
    % Dendrogram
    subplot('position',[0.0500    0.8178    0.9000    0.1322]);
    [H,~,perm] = dendrogram(ea_links{s,1},size(ea{s,1},1),'colorthreshold',th(s,1)); % dendrogram
    
    % Color Dendrogram By Cluster Colormap 
    lineColours = cell2mat(get(H,'Color')); % get line colours
    colourList = unique(lineColours,'rows'); % find unique line colours (starts with black)
    for c = 2:size(colourList,1) % for each colour (excluding black)
        i = ismember(lineColours, colourList(c,:),'rows'); % find lines of this color (logical)
        scrap = cell2mat(get(H(i),'XData')); % find the points that these lines cover
        scrap = scrap(:)'; % vectorise  
        scrap = unique(scrap(floor(scrap) == scrap)); % keep unique whole data points
        lineColours(i,:) = repmat(cmap_cluster{s,1}(mode(idx_numComp_sorted{s,1}...
            (sample_a{s,1}(perm(scrap)))),:),sum(i),1); % assign new color
    end
    for l = 1:size(H,1) % for each line
        set(H(l), 'Color', lineColours(l,:)); % set new color
    end
    axis off;
    
    % Maximum lifetime threshold line 
    line(get(gca,'xlim'), [th(s,1) th(s,1)],'LineStyle',':','LineWidth',1.5);
    text(1,double(th(s,1)),'Maximum Lifetime Cut','verticalalignment','bottom',...
        'FontName','Calibri','FontSize',18);
    
    % EA Matrix
    subplot('position', [0.0500    0.0500    0.9000    0.7478]);
    imagesc(ea{s,1}(perm,perm));
    axis off
    c = colorbar;
    c.Label.String = 'E.A. Index';
    c.Location = 'southoutside';
    c.FontSize = 18;
end

clear s H perm lineColours colourList c i l scrap

%% Bouts in PCA Space 
% figure; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
% set(gca,'FontName','Calibri'); hold on; 
%     
% for c = numComp(1):-1:1
%     scatter(X{1,1}(sample_a_n{1,1}(sample_a_n{1,1}(:,2) == c,1),1),...
%         X{1,1}(sample_a_n{1,1}(sample_a_n{1,1}(:,2) == c,1),2),'marker','.',...
%         'markerfacecolor',cmap_cluster{1,1}(c,:),...
%         'markeredgecolor',cmap_cluster{1,1}(c,:),...
%         'markerfacealpha',0.8,'markeredgealpha',0.8); 
% 
% end 

% Options: 
    % sample 1000 for each cluster and do some sort of density plot 
    % use sample_a_n points
    % tsne? 

%% Cluster Parameters - Fits 

tic
parameter_dists = cell(2,1); % structure {s,parameters}(cluster,min:max(parmater:))

for s = 1:2 % for active & inactive
    for p = 3:size(cells{s,1},2) % for each parameter
        for k = 1:numComp(s) % for each cluster
            clear pd;
            pd = fitdist(cells{s,1}(idx_numComp_sorted{s,1}==k,p),'kernel','Width',1); % Fit
            parameter_dists{s,p-2}(k,:) = pdf(pd,min(cells{s,1}(:,p)):max(cells{s,1}(:,p)));
        end
    end
end
toc

clear s p k pd 

%% Cluster Parameters - Fits Figure 

figure;
counter = 1; % start a counter 
for s = 1:2 % for active & inactive
    for p = 3:size(cells{s,1},2) % For each parameter
        % Figure Settings
        subplot(2,4,counter); hold on;
        box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
        if s ~= 2 % for active parameters 
            title(parameters{1}{p-2}); % Add title
        else % for inactive bout length 
            title(parameters{1}{10}); % Add title   
        end
        
        % Plot 
        crop = max(cells{s,1}(:,p)); 
        for k = 1:numComp(s) % for each cluster
            plot((min(cells{s,1}(:,p)):crop)/unit_conversion{1}(s,p-2),...
                parameter_dists{s,p-2}(k,:),'color',cmap_cluster{s,1}(k,:),'linewidth',3)
        end
        
        % Axes 
        try 
            set(gca,'XTick',...
                [min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2), (min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2))*10,...
                crop/unit_conversion{1}(s,p-2)]); % set x tick labels
        catch 
            set(gca,'XTick',...
                [min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2), (crop/unit_conversion{1}(s,p-2))/2,...
                crop/unit_conversion{1}(s,p-2)]); % set x tick labels
        end 
        axis([min(cells{s,1}(:,p))/unit_conversion{1}(s,p-2) crop/unit_conversion{1}(s,p-2) ...
            min(parameter_dists{s,p-2}(:)) max(parameter_dists{s,p-2}(:))]); % Set axis limits
        
        % Set decimal places depending on units
        if unit_conversion{1}(s,p-2) > 1
            xtickformat('%.2f');
        else
            xtickformat('%.0f');
        end
        
        set(gca,'XScale','log'); % set log axis
        xlabel(units{1}(p-2),'Fontsize',12); % X labels
        ylabel('Probability','Fontsize',12); % Y label
        
        counter = counter + 1; % add to counter 
        
        % Legend Plot 
        if counter == 8 
            subplot(2,4,counter); hold on;
            box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
            title('Legend'); 
            
            for s_2 = 1:2 % for active & inactive 
                for k = 1:numComp(s_2) % for each cluster 
                    scatter(k,mod(s_2,2),300,...
                        'markerfacecolor',cmap_cluster{s_2,1}(k,:),...
                        'markeredgecolor',cmap_cluster{s_2,1}(k,:));
                end 
            end 
            
            axis([0 max(numComp)+1 -1 2]);
            xlabel('Module Number','Fontsize',12); % X labels
            set(gca,'XTick',1:max(numComp)); 
            set(gca,'YTick',[0 1]); 
            set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12);
        end 
        
    end
end

clear counter s p crop k s_2

%% Bout Proportions 
    % 43 mins for 452 fish

tic
bout_proportions{1,1} = nan(max(fish_tags{1,1}),numComp(1),max(parameter_indicies{1,1}),...
    'single'); % fish x clusters x time windows 
bout_proportions{2,1} = nan(max(fish_tags{2,1}),numComp(2),max(parameter_indicies{2,1}),...
    'single'); % fish x clusters x time windows 

for s = 1:2 % for active & inactive  
    % note that comms overhead makes this faster as a for rather than a parfor loop  
    for f = 1:max(fish_tags{s,1}) % For each fish
        for c = 1:numComp(s) % For each bout type
            for t = 1:max(parameter_indicies{s,1}(fish_tags{s,1}==f)) % For each time window that fish uses 
                bout_proportions{s,1}(f,c,t) = sum(fish_tags{s,1}==f & idx_numComp_sorted{s,1}==c ...
                    & parameter_indicies{s,1}==t)/...
                    sum(fish_tags{s,1}==f & parameter_indicies{s,1}==t); 
                % the number of times fish (f) uses cluster (c) @ time (t) 
                % divided by the number of bouts fish (f) has @ time (t) 
                
                % Note - will return zero's when a fish doesn't use a
                % particular bout type :-)
            end
        end
    end
end
toc

clear s f c t 

%% Bout Proportion Stats

for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings 
    
    for s = 1:2 % for active & inactive
        
        % Grouping Variables
        anova_group = repmat(i_group_tags(i_experiment_reps==er),...
            [size([days{set_token} nights{set_token}],2),1])'; % groups
        anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
            [size([days{set_token} nights{set_token}],2),1])'; % experiments
        
        anova_time = [];
        for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
            anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
            % Allocate alternating zeros and ones to each time window
        end
        anova_time = anova_time';
        
        % Development Grouping Variable
        if size(days_crop{set_token}(days{set_token}),2) == ...
                size(nights_crop{set_token}(nights{set_token}),2) ...
                && size(days_crop{set_token}(days{set_token}),2) > 1 
            % If there are an equal number of windows 
            % & more than 1 time window 
            
            anova_development = []; % development
            anova_development = zeros(1,size(anova_group,2)); % Pre-allocate
            d = 1:size(anova_development,2)/(size(time_window{set_token}(1):...
                time_window{set_token}(2),2)/2):...
                size(anova_development,2); % divide into "24h" windows
            for t = 1:size(d,2)-1
                anova_development(d(t):d(t+1)-1) = t;
            end     
        else 
            anova_development = ones(size(anova_experiment)); % use all ones  
        end
        
        % Comparison
        for c = 1:numComp(s) % For each cluster
            clear scrap;
            scrap = permute(bout_proportions{s,1}(i_experiment_reps==er,...
                c,time_window{set_token}(1):time_window{set_token}(2)),[1 3 2]);
            scrap = scrap(:)'; % Vectorise
            
            [twa.bp.p{s,er}(:,c),~,twa.bp.stats{s,er,c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development,anova_experiment},...
                'display','off','model','full');
        end
        
        clear anova_development anova_experiment anova_group anova_time ... 
            scrap 
    end
end

clear er set_token s anova_group anova_experiment anova_time anova_development ...
    c scrap 