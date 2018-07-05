% Motif Examples Figure 

%% Options     
set(0,'DefaultFigureWindowStyle','docked'); % dock figures
set(0,'defaultfigurecolor',[1 1 1]); % white background

%% Module Distributions 

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'grammar_mat'); 

% calculate every motifs module distribution 
motif_m_dists = zeros(length(grammar_mat{1,1}),max(grammar_mat{1,1}(:))); % motifs x modules  
edges = 1:(max(grammar_mat{1,1}(:))+1); % hist count edges 
for s = 1:length(grammar_mat{1,1}) % for each motif 
    scrap = histcounts(grammar_mat{1,1}(s,:),edges); % count it's modules
    motif_m_dists(s,:) = scrap/sum(scrap); % normalise 
end 

% calculate each lengths module distribution (NaN for unseen lengths)
scrap = [sum(isnan(grammar_mat{1,1})==0,2) grammar_mat{1,1}]; % lengths & grammar 
for i = 1:max(scrap(:,1)) % for each motif length 
    toy = scrap(scrap(:,1) == i,2:end); % take motifs of this length   
    toy = toy(:); % vectorise 
    toy(isnan(toy)) = []; % remove nan values 
    temp = histcounts(toy,edges); % count modules  
    motif_a_dists(i,:) = temp/sum(temp); % normalise 
end 

clear edges s i toy temp

%% Representative Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'sleep_cells');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'idx_numComp_sorted');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'bouts');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'cmap_cluster_merge');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'numComp'); 

% mean inactive module length (frames) 
ibl = grpstats(sleep_cells(:,3),idx_numComp_sorted{2,1},'mean'); % average length  
ibl(1) = []; % remove NaN's 

%% Figure: Typical Motifs of Each Length V2 
figure; hold on;
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);

b = 1; % baseline counter (plots from bottom to top)
for i = flip(unique(scrap(:,1)))' % for each motif length (longest to shortest) 
    
    % Find a representative motif of this length
    temp = motif_m_dists;
    temp(scrap(:,1) ~= i,:) = NaN; % keep only motifs of this length
    
    if sum(isnan(temp(:,1)) == 0) >= 3 % if there are more than 3 motifs of this length
        idx = knnsearch(temp, motif_a_dists(i,:)); % find the motif that best
        % matches the average module stats for this length
        seq = grammar_mat{1,1}(idx,:); % find this motifs module sequence
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
        x_lims(b,:) = [1 a]; % track the x limits
        l_n(b) = i; % keep length 
        b = b + 1; % add to baseline counter
        
    end
    
end

x_lims = max(x_lims); 
axis([x_lims .5 (b-1)]); % hard coded axis 
set(gca,'XTick',round([1 x_lims(2)/2 x_lims(2)],2)); % hc x axis ticks  
set(gca,'XTickLabels',round([1 x_lims(2)/2 x_lims(2)]/25,2)); % hc x axis labels 
set(gca,'FontSize',32); 
xlabel('Time (Seconds)','Fontsize',32); 
ylabel({'Representative Motif' ; '(of Module Length)'},'Fontsize',32);
set(gca,'YTick',1:(b-1)); % hc y axis ticks 
set(gca,'YTickLabels',l_n,'Fontsize',32); % hc y axis labels  

%% Figure: Typical Motifs of Each Length V1 

figure; 

counter = 1; 
for i = unique(scrap(:,1))' % for each motif length
    ax(counter) = subplot(2,7,counter); hold on;
    set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);

    temp = motif_m_dists;
    temp(scrap(:,1) ~= i,:) = NaN; % keep only motifs of this length
    
    if sum(isnan(temp(:,1)) == 0) >= 3 % if there are more than 3 motifs of this length
    idx = knnsearch(temp, motif_a_dists(i,:)); % find the motif that best 
        % matches the average module stats for this length  
    seq = grammar_mat{1,1}(idx,:); % find this motifs module sequence 
    seq(isnan(seq)) = []; % remove nan values 
    a = 1; % start a counter (frames) 
    
    for t = 1:length(seq) % for each module in the sequence 
        if seq(t) <= numComp(1) % for the inactive modules 
            plot([a (a+ibl(seq(t)))],[0 0],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot 
            a = a + ibl(seq(t)); % add to time 
        else % for the active modules 
            plot(a:(a+length(nanmean(bouts{1,seq(t)-numComp(1)}))+1),...
                [0 nanmean(bouts{1,seq(t)-numComp(1)}) 0],...
                'color',cmap_cluster_merge(seq(t),:),'linewidth',5); % plot
            a = a + length(nanmean(bouts{1,seq(t)-numComp(1)})) + 1; % add to time 
        end
    end
    
    if counter ~= 1 % for most subplots turn the axis off 
        axis off;
    end
    
    x_lims(counter,:) = [1 a]; % track the x limits   
    y_lims(counter,:) = ylim; % track the y limits

    counter = counter + 1; % add to counter 
    
    end 
end

set(ax,'XLim',[-10 max(x_lims(:))]); % set x limits for all subplots    
set(ax,'Ylim',[0 27]); % set y limits for all subplots (hard coded)  

p = 1:(diff([1 max(x_lims(:))])/size(motif_a_dists,1)):diff([1 max(x_lims(:))]); 
    % calculate spacing for module overlay 
    
% Module Overlay 
counter = 1; % start a counter 
for i = unique(scrap(:,1))' % for each motif length
    
    temp = motif_m_dists;
    temp(scrap(:,1) ~= i,:) = NaN; % keep only motifs of this length
    if sum(isnan(temp(:,1)) == 0) >= 3 % if there are more than 3 motifs of this length
        
        idx = knnsearch(temp, motif_a_dists(i,:)); % find the motif that best
        % matches the average module stats for this length
        seq = grammar_mat{1,1}(idx,:); % find this motifs module sequence
        seq(isnan(seq)) = []; % remove nan values
        axes(ax(counter)); % move to appropriate axis 
        for s = 1:length(seq) % for each module
            scatter(p(s),26,90,...
                'markerfacecolor',cmap_cluster_merge(seq(s),:),...
                'markeredgecolor','k');
        end
        counter = counter + 1; % add to counter 
        
    end
end

% Axis 
%ax(1).XRuler.Axle.LineStyle = 'none'; % remove axis line 
%ax(1).YRuler.Axle.LineStyle = 'none'; % remove axis line 
set(ax(1),'XTick',max(x_lims(:)));  % set x limits 
set(ax(1),'XTickLabel',round(max(x_lims(:))/25,2)); % convert to seconds (2dp) 
set(ax(1),'YTick',25); % set y limits (hard coded)   
set(ax(1),'Fontsize',16); % font size 
ax(1).XLabel.String = 'Time (Seconds)'; % x label 
ax(1).YLabel.String = 'Delta Px'; % y label 

% "Legend"
axes(ax(1)); 

for s = 1:length(cmap_cluster_merge) % for each module 
    if s <= numComp(1) % for the inactive modules 
        scatter(max(x_lims(:))-40,s+16,90,'markerfacecolor',cmap_cluster_merge(s,:),...
            'markeredgecolor','k');
    else % for the active modules 
        scatter(max(x_lims(:)),s-numComp(1)+16,90,'markerfacecolor',cmap_cluster_merge(s,:),...
            'markeredgecolor','k');
    end
end 
text(max(x_lims(:))-65,23,'Module','FontName','Calibri','Fontsize',16);

%% Real Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'gCount_norm');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'threads'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_reps');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'i_experiment_tags')
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'states');
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'raw_data'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'fish_tags_cm'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'offset'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Draft_1\Post_Bout_Transitions.mat', 'uniqueSeqs'); 

counter = 1;
for i = unique(scrap(:,1))' % for each motif length
    subplot(2,7,counter); axis off
    if i ~= 11
        temp = motif_m_dists;
        temp(scrap(:,1) ~= i,:) = NaN;
        idx = knnsearch(temp, motif_a_dists(i,:));
        
        s = idx;
        % Plot Sequence - from randomly chosen fish
        % Randomly choose a fish who uses this motif
        exampleFish = datasample(find(sum(squeeze(gCount_norm{1,1}(idx,:,:))) >= 1),1);
        er = i_experiment_reps(exampleFish);
        Tlocs = datasample(strfind(threads{exampleFish,1,1}',uniqueSeqs{1,1}{s,1}),...
            1,'replace',false); % sample an instance of this sequence (position in threads)
        Rlocs = threads{exampleFish,2,1}(Tlocs,1):...
            threads{exampleFish,2,1}(Tlocs+(size(uniqueSeqs{1,1}{s,1},2)-1),2); % start-end (position in frames)
        Rlocs = Rlocs + offset(er,(i_experiment_tags(exampleFish) - ...
            (min(unique(i_experiment_tags(i_experiment_reps == er))) - 1))); % position in offset data
        if er ~= 1 % fish tag for raw_data
            exampleFish = exampleFish - fish_tags_cm(er - 1);
        end
        % ax = imagesc([1,size(raw_data{er,1}(exampleFish,Rlocs),2)],...
        %     [0,max(raw_data{er,1}(exampleFish,Rlocs))],...
        %     repmat(states{er,1}(exampleFish,Rlocs),[max(raw_data{er,1}(exampleFish,Rlocs)),1]));
        % hold on; colormap(cmap_cluster_merge); set(gca,'Ydir','Normal');
        plot(raw_data{er,1}(exampleFish,Rlocs),'k','linewidth',3); hold on;
        counter = counter + 1; 
    end
    
end
