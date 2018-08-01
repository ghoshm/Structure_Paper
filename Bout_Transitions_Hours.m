% Bout_Transitions_Hours

% marcus.ghosh.11@ucl.ac.uk 

%% Info 

%% Assumptions 

%% Dependencies 

% Files  
    % Compression Code 
        % https://github.com/aexbrown/Behavioural_Syntax/tree/master/compression
    
    % mi toolbox 
        % https://uk.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation?s_tid=prof_contriblnk
    
    % mrmr_d_matlab_src 
        % https://uk.mathworks.com/matlabcentral/fileexchange/14916-minimum-redundancy-maximum-relevance-feature-selection?s_tid=prof_contriblnk
          
% MATLAB Products 
    % Statistics and Machine Learning Toolbox 
    % Curve Fitting Toolbox 
    
%% Options     
set(0,'DefaultFigureWindowStyle','docked'); % dock figures
set(0,'defaultfigurecolor',[1 1 1]); % white background
 
%% Load data 
%load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Draft_1\180519.mat'); 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\Thesis\180725_2.mat'); % thesis 

%% Generate Hour Indexing Variable
    % Note that as this is it will work better for experiments where
    % you cut out a middle section

parameter_indicies_hours{1,1} = []; % allocate active
parameter_indicies_hours{2,1} = []; % allocate inactive

for e = 1:size(experiment_reps,2) % for each experiment
    
    gap = round((lb{e}(time_window{e}(2)+1) - lb{e}(time_window{e}(1)))/...
        (size(days{e},2)*24)); % average frames per hour 
    
    edges = [1 lb{e}(time_window{e}(1)):gap:...
        lb{e}(time_window{e}(2)+1) lb{e}(end,1)]; % hour times 
    
    if size(edges,2) ~= (size(days{e},2)*24) + 3 % if there are not the correct number of edges 
        edges = [1 lb{e}(time_window{e}(1)):gap:...
        lb{e}(time_window{e}(2)+1) lb{e}(time_window{e}(2)+1) lb{e}(end,1)]; 
        % add a bin with slighlty different number of frames
    end
    
    % Active
    parameter_indicies_hours{1,1} = [parameter_indicies_hours{1,1} ; ...
        single(discretize(wake_cells(experiment_tags{1,1} == e,1),...
        edges))]; % bin active data into hours 
    
    % Inactive
    parameter_indicies_hours{2,1} = [parameter_indicies_hours{2,1} ; ...
        single(discretize(sleep_cells(experiment_tags{2,1} == e,1),...
        edges))]; % bin inactive data into hours 
    
end

clear e gap edges

%% Threading Every Hour With Multiple Shuffles of Data 
% Very Slow, Roughly 7 hours for 452 fish, each with 10 shuffles 
% Note that I'm now keeping only the real time windows
% Rather than storing multiple copies of these (which are redundant)

% Shuffles
shuffles = 10; % hard coded number of shuffles

% Pre-allocation
threads = cell(max(fish_tags{1,1}),3,(1+shuffles)); % fish x (clusters,times (start & stop),time windows) x (samples & shuffled controls)
idx_numComp_sorted{1,1} = idx_numComp_sorted{1,1} + max(idx_numComp_sorted{2,1}); % Assign higher numbers to the wake bouts
idx_numComp_sorted{2,1}(isnan(idx_numComp_sorted{2,1})) = 0; % replace nan-values with zero

tic
for f = 1:max(fish_tags{1,1}) % For each fish
    % Pre-allocate
    % Data
    threads{f,1,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % clusters
    threads{f,2,1} = nan(size(threads{f,1,1},1),2,'single'); % times (start & stop)
    threads{f,3,1} = nan(size(threads{f,1,1}),'single'); % time windows
    
    % Deterime starting state (a = active, b = inactive)
    if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
        a = 1; b = 2; % start active
    else % If the fish starts inactive
        a = 2; b = 1; % start inactive
    end
    
    % Fill in Real Data
    % Clusters
    threads{f,1,1}(a:2:end,1) = idx_numComp_sorted{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active clusters
    threads{f,1,1}(b:2:end,1) = idx_numComp_sorted{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive clusters
    % Times
    threads{f,2,1}(a:2:end,1:2) = wake_cells...
        (fish_tags{1,1} == f,1:2); % Fill in active times
    threads{f,2,1}(b:2:end,1:2) = sleep_cells...
        (fish_tags{2,1} == f,1:2); % Fill in inactive times
    % Time Windows
    threads{f,3,1}(a:2:end,1) = parameter_indicies_hours{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active time windows
    threads{f,3,1}(b:2:end,1) = parameter_indicies_hours{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive time windows
    
    % Generate Shuffled Control Data
    % Note that each time window (e.g. hour 1 vs hour 2 etc) is
    % shuffled individually within each fish
    % Preserving the developmental & day/night cluster frequencies
    
    % Clusters
    for tc = 2:size(threads,3) % for each shuffle
        
        % Determine Starting State
        % Note that it's necessary to do this after each shuffle to
        % re-set the starting state
        if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
            a = 1; b = 2; % start active
        else % If the fish starts inactive
            a = 2; b = 1; % start inactive
        end
        
        % Shuffle data
        for tw = min(parameter_indicies_hours{1,1}(fish_tags{1,1} == f)):...
                max(parameter_indicies_hours{1,1}(fish_tags{1,1} == f)) % for each time window
            data = []; scrap = []; % empty structures
            
            % take real clusters from this time window
            scrap{1,1} = idx_numComp_sorted{1,1}(fish_tags{1,1} == f & parameter_indicies_hours{1,1} == tw,1); % active
            scrap{2,1} = idx_numComp_sorted{2,1}(fish_tags{2,1} == f & parameter_indicies_hours{2,1} == tw,1); % inactive
            
            % shuffle & merge the data 
            data = nan((size(scrap{1,1},1) + size(scrap{2,1},1)),1,'single'); % allocate (active + inactive clusters)
            data(a:2:end,1) = scrap{1,1}(randperm(length(scrap{1,1}))); % fill active clusters
            data(b:2:end,1) = scrap{2,1}(randperm(length(scrap{2,1}))); % fill inactive clusters
            
            threads{f,1,tc} = [threads{f,1,tc} ; data]; % store data
            
            % Deterime current final state (a = wake, b = sleep)
            if threads{f,1,tc}(end,1) <= numComp(2) % If the fish ends inactive
                a = 1; b = 2; % next will be active
            else % If the fish ends active
                a = 2; b = 1; % next will be inactive
            end
        end
        
    end
    
    % Crop out extra time windows
        % Introduces NaN values to threads 
    threads{f,3,1}(threads{f,3,1} == 1) = NaN; % remove first window 
        % Note that experiments whos first time window is the first frame
        % will have their first window assigned a value of 2
    threads{f,3,1}(threads{f,3,1} == max(threads{f,3,1})) = NaN; % remove last window 
    threads{f,3,1} = threads{f,3,1} - 1; % shift windows back one 
        
    % Report progress
    if mod(f,50) == 0 % every 50 fish
        disp(horzcat('Threaded Fish ',num2str(f),' of ',...
            num2str(max(fish_tags{1,1})))); % Report progress
    end
    
end
toc

clear f a b tc tw data scrap

save('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731_Hours.mat','-v7.3');

%% Load Data (180524 / 180731)

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Thesis\180731_Hours.mat'); 

%% Chunk Compression 
% Compress cluster chunks 
    % Note that threads{f,3,1} now contains hours and NaN values 
    % outside of the time windows of interest   

% settings
step = 500;
sMax = max(idx_numComp_sorted{1,1}) + 1; % Maximum states + 1 (first new symbol)
nMax = 10; % Maximum n-grams to consider

% allocate
chunks = cell(size(threads,1),1); % fish x 1
totSavings_cells = cell(size(threads,1),1); % fish x 1

tic
for f = 1:size(threads,1) % for each fish
    chunks{f,1} = find(threads{f,3,1} == 1,1,'first'):...
        step:find(isnan(threads{f,3,1}) == 0,1,'last');
    
    for t = 1:(length(chunks{f,1})-1) % for each chunk
        [~,~, totSavings_cells{f,1}(1,t)] = ...
            compressSequenceNFast(threads{f,1,1}(chunks{f,1}(1,t):(chunks{f,1}(t+1)-1))',...
            sMax,nMax); % compress returning totSavings
    end
    
    chunks{f,1}(end) = []; % remove the spare edge
    
    disp(num2str(f));
    
end
toc 

clear sMax nMax f t

%% Reshape Chunks & Calculate Compression  
    % The compressibility of a sequence of uncompressed length l is given by the sum of the savings S
    % at each iteration divided by l (Gomez-Marin et al.,2016)
    
% settings
tw = 48; % hard coded maximum number of time windows
dn_hour = ones(1,tw)*2; % day and night hours
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours

% allocate
compressibility = nan(size(threads,1),tw,'single'); % fish x max time hour bins

% fill data
for f = 1:size(threads,1) % for each fish
    for h = 1:max(threads{f,3,1}) % for each hour
        try
            compressibility(f,h) = nanmean(totSavings_cells{f,1}(...
                threads{f,3,1}(chunks{f,1}) == h)/step); 
            % take a mean score every hour  
        catch 
            % catches hours with less cluster than chunks 
            % note that these will remain NaN
        end
    end
    
end

clear f h 

%% Compressibility Every Hour - Figure
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings
figure;
hold on; set(gca,'FontName','Calibri'); clear scrap;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    clear data;
    data = compressibility(i_experiment_reps == er & i_group_tags == g,:);
    
    if er == 0 % for the WT experiments
        plot(data',...
            'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
        scrap(1,g) = min(data(:));
        scrap(2,g) = max(data(:));
        
    else 
        scrap(1,g) = min(nanmean(data) - nanstd(data));
        scrap(2,g) = max(nanmean(data) + nanstd(data)); 
    end
    
    legend_lines(g) = errorbar(nanmean(data),nanstd(data),...
        'color',cmap{set_token}(g,:),'linewidth',3);
    legend_cell{g} = horzcat(geno_list{set_token}.colheaders{g},', n = ',...
            num2str(size(data,1)));
end

y_lims = [(min(scrap(1,:)) - min(scrap(1,:))*0.05) ...
    (max(scrap(2,:)) + max(scrap(2,:))*0.05)]; % Add a bit of space either side
   
% Night Patches 
a = 1; night_start = 15; % hard coded counter
for n = 1:size(nights{set_token},2) % For each night
    r(a) = rectangle('Position',[(night_start) y_lims(1)...
        9 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    a = a + 1; night_start = night_start + 24; % Add to counters
end
    
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel('Time (Hours)','Fontsize',32); % X Labels 
ylabel({'Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels
axis([1 (n*24) y_lims]); 
if er == 0 % for the WT data
   [~,icons,plots,~] = legend(legend_lines,legend_cell,'Location','northwest');
else
   [~,icons,plots,~] = legend(legend_cell,'Location','northeast');
end
legend('boxoff'); 
set(icons(1:g),'Fontsize',32) ; set(plots,'LineWidth',3);

clear et set_token g data scrap legend_lines legend_cell y_lims a night_start n r icons plots   

%% Compressibility Day vs Night
    % Note 180612: make sure that dn_hour is set to 1-2
figure;
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    subplot(1,2,er); counter = 1; % counts groups for plots
    hold on; set(gca,'FontName','Calibri'); clear scrap;
    
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        clear data;
        % average day per fish
        data(:,1) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 1),2);
        % average night per night
        data(:,2) = nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,dn_hour == 2),2);
        plot([counter,counter+1],data,...
            'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
        errorbar([counter,counter+1],nanmean(data),nanstd(data),...
            'color',cmap{set_token}(g,:),'linewidth',3);
        counter = counter + 2; % add to counter
        
        scrap(1,g) = min(data(:));
        scrap(2,g) = max(data(:));
    end
    
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    if er == 0 % for the WT Data
        set(gca, 'XTick', [1 2]); % set X-ticks
        set(gca,'XTickLabels',{'Day','Night'}); % X Labels
    else % for the other experiments
        set(gca,'XTick',1.5:2:(max(i_group_tags(i_experiment_reps == er))*2)+.5);
        set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
    end
    ylabel({'Compressibility' ; '(per 500 modules)'},'Fontsize',32); % Y Labels
    axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 ...
        (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);
end

clear er set_token g scrap counter data

%% Compressibility Two Way ANOVA
dn_hour(1:14) = 1; dn_hour(15:24) = 2; dn_hour(25:38) = 3; dn_hour(39:48) = 4; 

for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    clear data; 
    
    for t = 1:max(dn_hour) % for each day/night 
        data(:,t) = nanmean(compressibility(i_experiment_reps == er,dn_hour == t),2); 
    end % average for each fish 
    data(:,find(sum(isnan(data)) == size(data,1),1,'first'):end) = []; % remove nans on the end 

    % Grouping Variables
    anova_group = repmat(i_group_tags(i_experiment_reps==er),...
        [size(data,2),1])'; % groups
    
    anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
        [size(data,2),1])'; % experiments
    
    anova_time = [];
    for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
        anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
        % Allocate alternating zeros and ones to each time window
    end
    anova_time = anova_time';
    
    % Development Grouping Variable
    if size(days_crop{set_token}(days{set_token}),2) == ...
            size(nights_crop{set_token}(nights{set_token}),2) ...
            && size(days_crop{set_token}(days{set_token}),2) > 1 % If there are an equal number of windows (>1)
        
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
    data = data(:)';
    
    [twa.rc.p{er}(:,1),~,twa.rc.stats{er}] = anovan(data,...
        {anova_group,anova_time,anova_development,anova_experiment},...
        'display','off','model','full');
    
end

clear er anova_group anova_experiment anova_time t anova_development d data

% Saved as (180525_Hours / 180801_Hours) 

%% -> Legion 

% load grammar freq data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Grammar_Hours_Results_Final.mat',...
    'gCount_norm');

% variables 
dn_hour = ones(1,size(gCount_norm{1,1},2))*2; % day and night hours 
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours 
shuffles = 11; 

%% Reformat Grammar Data 
temp = gCount_norm; clear gCount_norm; 

% Allocate 
for tc = 1:size(temp,2) % for real & shuffled data
    gCount_norm{1,tc} = zeros(size(temp{1,1},1),size(temp{1,1},2),size(temp,1),...
        'single'); % {real/shuffled} sequences x time windows x fish
end

% calculate inf replacements 
scrap = zeros(1,shuffles,'single'); % 1 x real & shuffled data  
scrap(1,1) = 1; 
inf_r = (scrap(1,1) - nanmean(scrap(1,2:end)))/nanstd(scrap); 

% Fill data 
for tc = 1:size(temp,2) % for real & shuffled data
    for f = 1:size(temp,1) % for each fish
        gCount_norm{1,tc}(:,:,f) = temp{f,tc}; % fill data
    end
    
    % Replace -Inf & Inf Values 
    gCount_norm{1,tc}(isinf(gCount_norm{1,tc}) & ...
        gCount_norm{1,tc} < 0) = -1*inf_r; % -ve inf 
    gCount_norm{1,tc}(isinf(gCount_norm{1,tc}) & ...
        gCount_norm{1,tc} > 0) = inf_r; % +ve inf 
      
end

clear temp tc scrap inf_r f shuffles

%% Load threaded data 

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180525_Hours.mat');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\180523.mat', 'grammar_mat');

%% Identifying Interesting Sequences (is)

% 180215 Notes
% 1. https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-mRMR_data.html#d119e2598
% 2. https://uk.mathworks.com/help/stats/fitcdiscr.html

% Settings
comps = 250; % number of sequences to identify

% mRMR Approach
for er = 1:max(experiment_reps) % for each experiment repeat
    
    % Grab Data
    if er == 1 % for the WT fish
        % Organised to be:
        % rows = fish 1 (all time windows), fish 2 (all time windows) etc 
        % columns = sequences 
        mRMR_data{er,1} = ...
            double(reshape(gCount_norm{1,1}(:,:,i_experiment_reps == er),...
            size(gCount_norm{1,1},1),[])');
        
        mRMR_tw{er,1} = repmat([1:size(gCount_norm{1,1},2)/2]',[sum(i_experiment_reps == er)*2,1]);
        
        mRMR_data{er,1}(mRMR_data{er,1} < 0) = 0; % Remove negative values for now
        
        mRMR_data{er,1}(mRMR_tw{er,1} > 14,:) = []; % keep only some data 
        mRMR_tw{er,1}(mRMR_tw{er,1} > 14,:) = []; % keep only some data 
        
    else
        % Organised to be:
        % rows = fish 1 (all time windows), fish 2 (all time windows) etc
        % columns = sequences
%         mRMR_data{er,1} = ...
%             double(reshape(gCount_norm{1,1}(:,time_window{set_token}(1):time_window{set_token}(2),...
%             i_experiment_reps == er),size(uniqueSeqs{1,1},1),[]))';
%         mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),1,...
%             size(time_window{set_token}(1):time_window{set_token}(2),2))';
%         mRMR_tw{er,1} = (mRMR_tw{er,1}(:)')';
    end
    
    % Pairwise Comparisons
    tic
    counter = 1; 
    for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group (hour)
        
        tags = ones(size(mRMR_tw{er,1}))*2; % all data 
        tags(mRMR_tw{er,1} == g_one) = 1; % group of interest 
        
        % mRMR
        [comps_v{er,1}(counter,:)] = mrmr_miq_d(...
            zscore(mRMR_data{er,1}),tags,comps);

        % Classifiers
        for s = 1:comps % for each comp sequence
            % Fit a linear classifier as you add features
            % Using 10 fold cross validation
            % Hold 10% of mRMR_data back by default
            Mdl = fitcdiscr(...
                zscore(mRMR_data{er,1}(:,...
                comps_v{er,1}(counter,1:s))),...
                tags,...
                'DiscrimType','linear','CrossVal','on');
            Mdl_loss{er,1}(counter,s) = kfoldLoss(Mdl);
            Mdl_loss{er,2}(counter,s) = nanstd(kfoldLoss(Mdl,'Mode','individual')); 
        end
        
        % Minimal Feature Space
        mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
            min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
                
        disp(num2str(counter)); 
        counter = counter + 1; 
    end
    disp(horzcat('Finished mRMR Comparisons ',num2str(er),' of ',...
        num2str(max(experiment_reps)))); % report progress
    toc
    
end

clear er set_token s Mdl

%% Interesting Motifs Figure
% Settings
er = 1; % set interest
set_token =  find(experiment_reps == er,1,'first'); % settings

figure; hold on;

% Motifs
clear scrap;
scrap = grammar_mat{1,1}(comps_v{er,1}(:,1),:); % grab motifs
scrap(:,sum(isnan(scrap)) == size(comps_v{er,1},1)) = [];
ax = imagesc(scrap,'AlphaData',isnan(scrap)==0); % imagesc with nan values in white
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(ax,'CDataMapping','direct');
c = colorbar; c.Label.String = 'Cluster';
xlabel('Position in Motif','Fontsize',32);
ylabel('Motif','Fontsize',32);

clear er set_token scrap ax c

%% Minimal Feature Space tSNE 
er = 1; % set interest
set_token = find(experiment_reps == er,1,'first'); % settings

if er == 1
    motifs = []; 
    for m = 1:size(mRMR_ms,2) % for each comparison
        motifs = [motifs comps_v{er,1}(m,1:mRMR_ms(er,m))];
    end 
    motifs = unique(motifs); 
  
    mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,motifs),...
        'Algorithm','barneshut','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
        'Perplexity',30,'Standardize',1,'Verbose',0);
else
%     scrap = [];
%     for c = 1:size(comps_v{er,1},1) % for each comparison
%         scrap = [scrap comps_v{er,1}(c,1:mRMR_ms(er,c))]; 
%     end
%     scrap = unique(scrap); 
%     mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,scrap),...
%         'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
%         'Perplexity',30,'Standardize',1,'Verbose',0);
end

% colormap
n = length(min(mRMR_tw{er,1}):max(mRMR_tw{er,1})); % number of colours 
CT = [linspace(cmap_2{1,1}(1,1),cmap_2{1,1}(2,1),n)'...
    linspace(cmap_2{1,1}(1,2),cmap_2{1,1}(2,2),n)'...
    linspace(cmap_2{1,1}(1,3),cmap_2{1,1}(2,3),n)']; 

figure; hold on;
for g = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
    if er == 1 % for the wt mRMR_data
        scatter(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1),mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2),...
            'markerfacecolor',CT(g,:),...
            'markeredgecolor',CT(g,:));
        %pause(3);
    else
        scatter(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1),mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2),...
            'markerfacecolor',cmap{set_token}(g,:),...
            'markeredgecolor',cmap{set_token}(g,:));
    end
    
end

%% Motifs Over Time 

% Best motif for each hour 
scrap = gCount_norm{1,1}(comps_v{er,1}(:,1),:,i_experiment_reps == er);

% Figure 
figure; hold on; 
for g = 1:max(mRMR_tw{er,1}) % for each group 
    plot(nanmean(scrap(g,:,:),3),'color',CT(g,:));
    pause(3);
end

%% Motifs By Hour Figure 
er = 1; % set interest
set_token = find(experiment_reps == er,1,'first'); % settings

figure; hold on;
for g = 1:max(mRMR_tw{er,1}) % for each group
    data = [squeeze(scrap(g,1:24,:))' ; squeeze(scrap(g,25:end,:))' ];
%     errorbar(nanmean(data),nanstd(data)/sqrt(124),...
%         'color',CT(g,:),'linewidth',3);
plot((nanmean(data) - min(nanmean(data)))./range(nanmean(data)),'color',CT(g,:),'linewidth',3)

    pause(3);
end

%% Variables for Motif Plotting 

ibl = grpstats(sleep_cells(:,3),idx_numComp_sorted{2,1},'mean');
ibl(1) = []; % remove NaN's
ibl(end) = 250; % crop longest module to 10s (hard coded)

cmap_cluster_merge = [cmap_cluster{2,1} ; cmap_cluster{1,1}]; 

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'bouts');

%% Representative Temporal Motifs 

er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings
figure; 

% Motifs 
subplot(1,2,1); 
hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Fontsize',32); 

b = 1; % baseline counter (plots from bottom to top)
for s = [comps_v{er,1}(1,1) 5933 26911] % for each motif (HARD CODED)
    seq = grammar_mat{1,1}(s,:); % find this motifs module sequence
    seq(isnan(seq)) = []; % remove nan values
    a = 1; % start a counter (frames)

    for t = 1:length(seq) % for each module in the sequence
        if seq(t) <= numComp(1) % for the inactive modules
            plot([a (a+ibl(seq(t)))],[b b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + ibl(seq(t)); % add to time
        else % for the active modules
            plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
                [b ((nanmean(bouts{1,seq(t)-numComp(1)})/28)+b) b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time
        end
    end
    
    x_lims(b) = a; 
    
    b = b + 1; % add to baseline counter
end

x_lims = max(x_lims); 
x_lims = ceil(x_lims/25)*25; 
axis([1 x_lims .75 b]); % hard coded axis
set(gca,'Layer','top');
set(gca,'XTick',[1 x_lims/2 x_lims]); % hc x axis ticks
set(gca,'XTickLabels',[1/25 (x_lims/2)/25 x_lims/25]); % hc x axis labels
xlabel('Time (Seconds)','Fontsize',32);
ylabel('Example Motif','Fontsize',32);
set(gca,'YTick',1:3);
set(gca,'YTickLabels',{'Startle','Night','Day'},'Fontsize',32); 

% Traces 
subplot(1,2,2); 
hold on; set(gca,'FontName','Calibri'); box off; set(gca,'Fontsize',32); 
s = comps_v{er,1}(1,1); % startle motif  
plot(rescale(nanmean(squeeze(gCount_norm{1,1}(s,:,i_experiment_reps == er)),2))+1,...
    'color',[1 .5 0],'linewidth',3); 
s = 5933; % night motif  
plot(rescale(nanmean(squeeze(gCount_norm{1,1}(s,:,i_experiment_reps == er)),2))+2,...
    'color',cmap_2{set_token}(2,:),'linewidth',3)
s = 26911; % day motif 
plot(rescale(nanmean(squeeze(gCount_norm{1,1}(s,:,i_experiment_reps == er)),2))+3,...
    'color',cmap_2{set_token}(1,:),'linewidth',3)

y_lims = [.5 b+.5];

% Night Patches
a = 1; night_start = 15; % hard coded counter
for n = 1:2
    r(a) = rectangle('Position',[(night_start) y_lims(1)...
        9 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    night_start = 39;  
end

axis([1 size(gCount_norm{1,1},2) .75 b]); % hard coded axis
set(gca,'Layer','top');
xlabel('Time (Hours)','Fontsize',32);
set(gca,'YTick',[]); 
ylabel('Z-Score (Rescaled)','Fontsize',32); 

%% Startle Motifs (Multiple)
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings

s = comps_v{er,1}(1,1); 
idx = knnsearch(nanmean(gCount_norm{1,1}(:,:,i_experiment_reps == er),3),...
    nanmean(gCount_norm{1,1}(s,:,i_experiment_reps == er),3),...
    'K',10);

scrap = gCount_norm{1,1}(idx,:,i_experiment_reps == er);
scrap_cmap = lbmap(length(idx),'BlueGray'); % colormap 

%% Startle Figure 
figure; hold on; set(gca,'FontName','Calibri'); box off; 
for s = 1:length(idx)
    data = squeeze(scrap(s,:,:))'; % fish x hours  
    errorbar(nanmean(data),nanstd(data)/sqrt(size(scrap,3)),...
        'color',scrap_cmap(s,:),'linewidth',3);
end

y_lims = ylim;

% Night Patches
a = 1; night_start = 15; % hard coded counter
for n = 1:2
    r(a) = rectangle('Position',[(night_start) y_lims(1)...
        9 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    night_start = 39;  
end

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel('Time (Hours)','Fontsize',32); % X Labels 
ylabel('Z-Score','Fontsize',32); % Y Labels
axis([1 (n*24) y_lims]); 

% Insert 
ax1 = axes('Position',[0.15 0.5 0.1 0.4]); hold on; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16); set(gca,'FontName','Calibri'); % Set Font
set(gca,'color','none'); 

% Motif
b = 1; % baseline counter (plots from bottom to top)
for s = flip(comps_v{er,1}(:,1))' % for each motif
    seq = grammar_mat{1,1}(s,:); % find this motifs module sequence
    seq(isnan(seq)) = []; % remove nan values
    a = 1; % start a counter (frames)

    for t = 1:length(seq) % for each module in the sequence
        if seq(t) <= numComp(1) % for the inactive modules
            plot([a (a+ibl(seq(t)))],[b b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + ibl(seq(t)); % add to time
        else % for the active modules
            plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
                [b ((nanmean(bouts{1,seq(t)-numComp(1)})/28)+b) b],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time
        end
    end

    b = b + 1; % add to baseline counter
end

clear er set_token s idx scrap scrap_cmap data y_lims a night_start n r ax1 ax    