%% Bout_Clustering Fit 

%% Preparing Model 

load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat', 'wake_cells');
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat', 'fish_tags');

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

mu = mean(X{1,1}); % mean Model feature values 

clearvars -except mu

%% Load Clustered Data  
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\Pre.mat','wake_cells','sleep_cells'); 
active = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_1.mat','idx');
inactive = load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Draft_1\180515_2.mat','idx'); 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'numComp')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'sleep_cells_nan_track')

% merge variables 
idx{1,1} = active.idx; idx{2,1} = inactive.idx; 

clear active inactive

% Remove NaN Values from Sleep matracies  
sleep_cells(sleep_cells_nan_track,3) = NaN; % sleep cells 
idx{2,1}(sleep_cells_nan_track,1) = NaN; % cluster assignments  

% Merge Variables for ease of Looping 
cells{1,1} = wake_cells; 
cells{2,1} = sleep_cells; % cells 

%% Sorting Clusters by Mean Length 
             
for s = 1:2 % for active & inactive
    mean_cluster_length = nan(1,numComp(s),'single'); % pre-allocate
    for c = 1:numComp(s) % For each cluster
        mean_cluster_length(c) = nanmean(cells{s,1}(idx{s,1}==c,3));
        % Calculate mean bout length
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    cluster_order(s,:) = O; 
    
    clear s mean_cluster_length c O;
end

clearvars -except mu cluster_order  

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
experiment_reps = [1 1 2 2 2]; % experiment groupings 
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

%% Preparing New Data   

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

% Load Model 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'coeff');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'score');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'knee_dim'); 

% Fit To Model's PCA Basis 
% https://stackoverflow.com/questions/13303300/how-to-project-a-new-point-to-pca-new-basis#
X{1,1} = bsxfun(@minus,X{1,1}, mu); % subtract old means from new data 
X{1,1} = X{1,1}*coeff; % project new data into pca basis 
X{1,1} = X{1,1}(:,1:knee_dim); % crop to number of PC's 

% Inactive 
sleep_cells_nan_track = isnan(sleep_cells(:,3)); % store nan locations  

X{2,1} = []; % empty X
X{2,1} = sleep_cells(:,3); 

% Convert from frames to seconds  
X{2,1}(experiment_tags{2,1} <= find(experiment_reps == 1,1,'last'),1) = ... 
    X{2,1}(experiment_tags{2,1} <= find(experiment_reps == 1,1,'last'),1)/fps{1}; 

X{2,1}(experiment_tags{2,1} > find(experiment_reps == 1,1,'last'),1) = ... 
    X{2,1}(experiment_tags{2,1} > find(experiment_reps == 1,1,'last'),1)/fps{end}; 

%% Assigning New Data to Clusters 

% Load variables 
nn = 50; 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'sample_a_n');
Mdl = load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'X'); 
Mdl = Mdl.X; % old data in PCA / Frame space 
Mdl{2,1} = Mdl{2,1}/25; % convert from Frames to seconds  

for s = 1:2 % for active and inactive
    disp('Fitting data into clusters'); % report progress
    
    ea_idx_n = sample_a_n{s,1}(:,2); % normalised cluster indicies
    
    idx{s,1} = zeros(size(X{s,1},1),1,'single'); % allocate
    cks = [1:round(size(X{s,1},1)/1000):size(X{s,1},1) (size(X{s,1},1)+1)]; % break into 1000 chunks
    for i = 1:(length(cks) - 1) % for each chunk
        kn = knnsearch(Mdl{s,1}(sample_a_n{s,1}(:,1),:),X{s,1}(cks(i):(cks(i+1)-1),:),...
            'K',nn); % find nn nearest neighbours
        idx{s,1}(cks(i):(cks(i+1)-1),1) = ...
            mode(reshape(ea_idx_n(kn(:)),size(kn)),2); % assign to mode neighbour cluster
    end
end

% "Repair Data" 
idx{2,1}(sleep_cells_nan_track == 1,1) = NaN; % assign NaN Values

% Convert back to frames
X{2,1}(experiment_tags{2,1} <= find(experiment_reps == 1,1,'last'),1) = ... 
    X{2,1}(experiment_tags{2,1} <= find(experiment_reps == 1,1,'last'),1)*fps{1}; 

X{2,1}(experiment_tags{2,1} > find(experiment_reps == 1,1,'last'),1) = ... 
    X{2,1}(experiment_tags{2,1} > find(experiment_reps == 1,1,'last'),1)*fps{end}; 

Mdl{2,1} = Mdl{2,1}*25; % convert from Frames to seconds  

%% Sorting Clusters by Model Length 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat', 'numComp'); 

for s = 1:2 % for active & inactive
    
    O = cluster_order(s,:);
    
    idx_numComp_sorted{s,1} = nan(size(idx{s,1},1),1,'single'); % Pre-allocate
    
    for c = 1:numComp(s)  % For each cluster
        idx_numComp_sorted{s,1}(idx{s,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    clear s O c
end

clear idx

% Data Saved Here (170725)

%% Bout Proportions 
    % ~= 70 mins for 443 fish

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
