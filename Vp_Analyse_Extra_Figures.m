%% Vp_Analyse Extra Figures

%% Figure - Activity (Cropped) 
    % Crops the activity so that only points with data from all experiments
        % are shown 

figure; hold on; clear legend_lines; set(gca,'FontName','Calibri'); % Set Font  
for e = 1:max(experiment_tags) % For each experiment
    for g = 1:max(group_tags) % For each group
        
        legend_lines(g) = shadedErrorBar(lb_sec(time_window(1)):...
            (lb_sec(time_window(2)+1)-1),nanmean(delta_px_sq_sec_smooth...
            (lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),group_tags == g & experiment_tags == e),2),...
            nanstd(delta_px_sq_sec_smooth(lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),...
            group_tags == g & experiment_tags == e)')/sqrt(size(find(group_tags == g & experiment_tags == e),1)),...
            'lineprops',{'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5))});
        
        legend_cols(g) = legend_lines(g).mainLine; % Store color
        
        legend_cell{g} = horzcat(geno_list.colheaders{g},', n = ',...
            num2str(size(find(group_tags == g),1)));
        % Append the group size to each group name
        
        if g == max(group_tags) && e == max(experiment_tags) % After the last group
            a = 1; % Start counter
            for n = 1:size(nights,2) % For each night
                y_lims = ylim; % Find the axis limits
                r(a) = rectangle('Position',[lb_sec(nights_crop(nights(n))) 0,...
                    (lb_sec(nights_crop(nights(n))+1)-1) - lb_sec(nights_crop(nights(n))) y_lims(2)],...
                    'FaceColor',night_color,'Edgecolor',[1 1 1]);
                uistack(r(a),'bottom'); % Send to back
                a = a + 1; % Add to counter
            end
        end
    end
end 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
[~,icons,plots,~] = legend(legend_cols,legend_cell,'Location','northwest');
legend('boxoff'); 
set(icons(1:max(group_tags)),'Fontsize',32) ; set(plots,'LineWidth',3);
axis([find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'first') + time_bins...
    find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'last')- time_bins 0 y_lims(2)]);

x = find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'first') + time_bins:...
    (60*60*12):find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'last')- time_bins; 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12})
xlabel('Time (Hours)','Fontsize',32);
ylabel('Delta Px','Fontsize',32);

clear a g h icons plots str legend_cell legend_cols legend_lines n r x y_lims 


%% Example Parameters 

figure;

% Number of Bouts 
subplot(1,2,1); hold on; col = 1; counter = 1;
clear data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font

p = 7;
col = 1;
for g = 1:max(group_tags) % for each group
    clear data;
    % Day
    data{1} = parameter_comparisons{p}(:,g,days_crop(days)); % extract data
    data{1}(isnan(data{1})) = []; % remove nan values
    data{1} = nanmean(reshape(data{1},[group_sizes(g),size(days_crop(days),2)]),2)/...
        unit_conversion(1,p); % take a mean & convert units
    
    % Night
    data{2} = parameter_comparisons{p}(:,g,nights_crop(nights));
    data{2}(isnan(data{2})) = [];
    data{2} = nanmean(reshape(data{2},[group_sizes(g),size(nights_crop(nights),2)]),2)/...
        unit_conversion(1,p);
    
    % Plot
    for e = 1:max(experiment_tags)
        spread_cols = plotSpread(data{1}(experiment_tags(group_tags == g) == e),...
            'xValues',g,'distributionColors',...
            cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),'showMM',2);
        spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
        spread_cols = plotSpread(data{2}(experiment_tags(group_tags == g) == e),...
            'xValues',g+1,'distributionColors',...
            cmap_2(col+1,:)+(1-cmap_2(col+1,:))*(1-(1/e^.5)),'showMM',2);
        spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
    end
    set(findall(gca,'type','line'),'markersize',30); % change marker sizes
    
    col = col + 2;
    
end
ylabel('Number of Bouts','Fontsize',32); % Y labels
set(gca,'xtick',1:2); % Set x ticks
set(gca,'xticklabel',{'Day','Night'},'Fontsize',32); % Name each group
xlabel('Time','Fontsize',32); % X labels

%% Inactive Bout Length (Binned) 
p = 10; 

subplot(1,2,2); hold on; col = 1; counter = 1;
clear scrap data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font

% V3
crop = 750; % crop

% Take day & night distributions across time windows
scrap{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
    [],size(parameter_dists{p}(:,:,days_crop(days)),2));
scrap{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
    [],size(parameter_dists{p}(:,:,nights_crop(nights)),2));
data{1,1} = scrap{1,1}(:,1:crop); 
data{1,2} = scrap{1,2}(:,1:crop); 
data{1,1}(:,end+1) = sum(scrap{1,1}(:,(crop+1):end),2); 
data{1,2}(:,end+1) = sum(scrap{1,2}(:,(crop+1):end),2); 

% Expand group tags to account for multiple days/nights
data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]);
data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]);

% Expand experiment tags to account for multiple days/nights
data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);

for e = 1:max(experiment_tags) % For each experimet
    col = 1;
    for g = 1:max(group_tags) % For each group
        for t = 1:2 % For day/night
            
            legend_lines(col) = shadedErrorBar((min(dist_boundaries(:,p)):crop)/unit_conversion(1,p),...
                nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,1:end-1)),...
                nanstd(data{1,t}(data{2,t} == g & data{3,t} == e,1:end-1))...
                ,'lineprops',{'color',cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5))});
            % Plot - note that this is now in appropriate units
            % (eg.seconds)
            
            legend_cols(col) = legend_lines(col).mainLine; % Save color     
            
    % Append time window to group names
            if e == 1 % For the first experiment
                if first_night == 2 % If starting in the day
                    if t == 1 % For the days
                        legend_cell{col} = horzcat(geno_list.colheaders{g},' - Day'); % Tag
                    else % For the nights
                        legend_cell{col} = horzcat(geno_list.colheaders{g},' - Night'); % Tag
                    end
                else % If starting in the night
                    if t == 1 % For the nights
                        legend_cell{col} = horzcat(geno_list.colheaders{g},' - Night'); % Tag
                    else % For the days
                        legend_cell{col} = horzcat(geno_list.colheaders{g},' - Day'); % Tag
                    end
                end
            end
            % Determine Axis boundaries
            y_lims(1,counter) = max(ylim);
            y_lims(2,counter) = 0;
            
            col = col + 1; % Add to color
            counter = counter + 1; % Add to counter
        end
        
    end
end
axis([min(dist_boundaries(:,p))/unit_conversion(1,p) crop/unit_conversion(1,p) ...
    min(y_lims(2,:)) max(y_lims(1,:))]); % Set axis limits
try
    set(gca,'XTick',...
        [min(dist_boundaries(:,p))/unit_conversion(1,p), (min(dist_boundaries(:,p))/unit_conversion(1,p))*10,...
        crop/unit_conversion(1,p)]); % set x tick labels
catch
    set(gca,'XTick',...
        [min(dist_boundaries(:,p))/unit_conversion(1,p),(0.5*crop)/unit_conversion(1,p),...
        crop/unit_conversion(1,p)]); % set x tick labels
end
% Set decimal places depending on units
if unit_conversion(1,p) > 1
    xtickformat('%.2f');
else
    xtickformat('%.0f');
end

set(gca,'XScale','log'); % set log axis
xlabel('Inactive Bout Length (Seconds)','Fontsize',32); % X labels
ylabel('Probability','Fontsize',32); % Y label

[~,icons,~,~] = legend(legend_cols,legend_cell,'Location','southeast');
legend('boxoff');
set(icons(1:max(group_tags)*2),'Fontsize',16);
set(icons((max(group_tags)*2)+1:2:end),'LineWidth',3);

%% Insert 
ax1 = axes('Position',[0.8 0.5 0.1 0.4]); hold on; 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font

col = 1; 
for t = 1:2 % for day/night 
    spread_cols = plotSpread(data{1,t}(experiment_tags(group_tags == g) == e,end),...
        'xValues',crop/unit_conversion(1,p),'distributionColors',...
        cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),'showMM',2);
    spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
    set(findall(gca,'type','line'),'markersize',15); % change
    col = col + 1; 
end 

% remember to add in > sign & change the line width of the SEB mean to 3

clear col g p spread_cols t xl icons data
