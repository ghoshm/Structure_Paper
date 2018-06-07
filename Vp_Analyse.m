% Vp_Analyse 

% marcus.ghosh.11@ucl.ac.uk 

%% Info 

% Analyses data output from Vp_Extract.m 

% Input: 
    % Workspace or workspaces output from Vp_Extract.m  

% Output: 
    % Figures & Statistics 
    
%% Assumptions 
 
% Experiments to be merged are 'similar': 
    % Groups are numbered the same across experiments (e.g. WT,Het,Hom - 1,2,3) 
    % Experiments have the same number of days/nights 
    % Experiments have the same groups? 

%% Dependencies 

% Files 
    % PlotSpread
        %https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points--beeswarm-plot-

    % ShadedErrorBar
        %https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
 
% MATLAB Products 
    % Neural Network Toolbox 
    % Image Processing Toolbox 
    % Statistics and Machine Learning Toolbox
    % Curve Fitting Toolbox 
    
%% Options
    
et = 1; % Analyse repeats as one experiment (1) or treat repeats seperatly (2) 
ct{1} = [2 3]; % days of interest 
ct{2} = [2 3]; % nights of interest 
ac = 1; % adjust color scheme (1) or not (2)  
set(0,'DefaultFigureWindowStyle','docked'); % dock figures
set(0,'defaultfigurecolor',[1 1 1]); % white background

%% Selecting Files  
 
%Load Data - using multiselect
 [filename, pathname] = uigetfile('*.mat', 'Select file or files','MultiSelect','on'); %Select files
 if isequal(filename,0) %If no file is selected
     error('No Files Selected') %Show Error
 else %If selected
     disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
 end
 
%% Combining Experiments (< Parameter Extracted Data)

 % Convert filename to a cell (useful for single files) 
 filename = cellstr(filename); 
 
 experiment_tags = [];
 parameter_time_sec_smooth = [];
 lb_frames = []; 
 lb_sec = [];
 group_tags = [];
 wake_cells = [];
 sleep_cells = [];
 top_up_bin = []; 
 parameter_indicies = [];
 parameter_matrix = [];
 parameter_comparisons = cell(1,12);
 for f = 1:size(filename,2) %For each file
     clear experiment; 
     experiment = load(strcat(pathname,filename{f})); %Load the mat file
     experiment_tags = [experiment_tags ; ones(size(experiment.group_tags,1),1)*f]; 
     group_tags = [group_tags ; experiment.group_tags]; 
             % Allocate experiment tags 
     
     % Merge variables 
     parameter_time_sec_smooth = [parameter_time_sec_smooth ; experiment.parameter_time_sec_smooth];  
     lb_sec = [lb_sec experiment.lb_sec]; 
     lb_frames = [lb_frames experiment.lb]; 
     wake_cells = [wake_cells experiment.wake_cells]; 
     sleep_cells = [sleep_cells experiment.sleep_cells]; 
     parameter_indicies = [parameter_indicies experiment.parameter_indicies]; 
     parameter_matrix = [parameter_matrix ; experiment.parameter_matrix]; 
     delta_px_sq_sec_smooth{f,1} = experiment.delta_px_sq_sec_smooth; 
     try
         top_up_bin{1,f} = experiment.top_up_bin; 
     catch
     end
     
     % Nab variables 
     if f == 1 % For the first file 
         parameters = experiment.parameters; 
         time_window = experiment.time_window; 
         time_bins = experiment.time_bins; 
         unit_conversion = experiment.unit_conversion; 
         cmap = experiment.cmap; 
         cmap_2 = experiment.cmap_2; 
         night_color = experiment.night_color; 
         nights = experiment.nights; 
         nights_crop = experiment.nights_crop; 
         geno_list = experiment.geno_list; 
         units = experiment.units; 
         units_2 = experiment.units_2; 
         parameter_smooth = experiment.parameter_smooth; 
         first_night = experiment.first_night; 
         days_crop = experiment.days_crop; 
         days = experiment.days;
     end 
        
     for p = 1:size(parameters,2) % For each parameter 
        parameter_comparisons{p} = [parameter_comparisons{p} ; ...
            experiment.parameter_comparisons{p}]; 
     end 
     
 end
 
 clear pathname filename geno f;

%% Aligning data

% Calculate time off-sets between repeats of experiments 
[~,c] = find(lb_sec == max(lb_sec(2,:)));
    % Find the column (c) with the longest starting window (in seconds) 
for e = 1:size(lb_sec,2) % For each experiment
    offset(1,e) = lb_sec(2,c) - lb_sec(2,e); % Calculate the offset
end

% Parameter time
parameter_time_sec_smooth(end+1,:) = cell(1,size(parameters,2)); 
start = max(lb_sec(2,:));
for p = 1:size(parameters,2) % For each parameter 
    parameter_time_sec_smooth{end,p} = nan(max(lb_sec(end,:) + offset),...
        size(group_tags,1),'single'); % Pre-allocate  
    
    for f = 1:max(experiment_tags) % For each experiment 
        parameter_time_sec_smooth{end,p}(start-lb_sec(2,f)+1:...
            start-lb_sec(2,f)+lb_sec(end,f),experiment_tags==f) = ...
            parameter_time_sec_smooth{f,p}; 
    end 
end 
parameter_time_sec_smooth = parameter_time_sec_smooth(end,:); % Keep merged data 

% Activity 
delta_px_sq_sec_smooth{end+1,1} = nan(max(lb_sec(end,:) + offset),...
    size(group_tags,1),'single'); % Pre-allocate
for f = 1:max(experiment_tags) % For each experiment 
    delta_px_sq_sec_smooth{end,1}(start-lb_sec(2,f)+1:...
        start-lb_sec(2,f)+lb_sec(end,f),experiment_tags==f) = ...
        delta_px_sq_sec_smooth{f,1}; 
end 
delta_px_sq_sec_smooth = delta_px_sq_sec_smooth{end,1}; % Keep merged data 

clear offset 

%% Normalise for time (% of time active)
lb_frames_diff = diff(lb_frames); % diff to get frames per window
top_up_counter = 1; % start counter

for g = 1:max(group_tags) % for each group
    
    clear list; counter = 1; 
    list = find(isnan(parameter_comparisons{1,8}(:,g,1))==0); % find values 

    for f = find(group_tags == g)' % for each fish in this group 
        for t = 1:size(parameter_comparisons{1,8},3) % for each time window
            if sum(isnan(sleep_cells{1,f}(parameter_indicies{2,f} == t,3))) == 0 % if there are no NaN Values
                parameter_comparisons{1,8}(list(counter),g,t) = ...
                    (parameter_comparisons{1,8}(list(counter),g,t)/...
                    lb_frames_diff(t,experiment_tags(f)))*100;
            else % If there are NaN values (water topped up)
                % Subtract the time taken to top the water up
                parameter_comparisons{1,8}(list(counter),g,t) = ...
                    (parameter_comparisons{1,8}(list(counter),g,t)/...
                    (lb_frames_diff(t,experiment_tags(f)) - ...
                    diff(top_up_bin{1,experiment_tags(f)}(top_up_counter,:))))*100;
                top_up_counter = top_up_counter + 1; % add to counter
            end
        end
        counter = counter + 1; % add to counter
        top_up_counter = 1; % re-set counter
    end
end
clear g f t counter top_up_counter list 

unit_conversion(1,8) = 1; units{8} = '% of Time'; % swap hours to % 
unit_conversion(1,12) = 1; units{12} = '% of Time'; % swap hours to % 

%% Other 

% Store time windows 
lb_sec = lb_sec(:,c); % Keep longest time windows 

% Re-determine Group sizes
for g = 1:max(group_tags) % For each group
    group_sizes(g) = size(find(group_tags == g),1);
end

% Handling Experiment Tags 
if et == 1
    experiment_tags(:,1) = 1;
    % If just using a single experiment
    %experiment_tags = ones(size(group_tags));
end

% Adjust color scheme 
if ac == 1
    if max(group_tags) == 1 % for WT experiments
        cmap(1,:) = [135 206 250]/255; % light sky blue
        cmap_2(1,:) = cmap;
        cmap_2(2,:) = [25 25 112]/255; % midnight blue
    else
        cmap_2 = flip(cmap_2); % flip cmap
        for c = 1:2:size(cmap_2,1)
            cmap_2([c c+1],:) = cmap_2([c+1 c],:); % swap the colors around
        end
        cmap = cmap_2(1:2:size(cmap_2,1),:); % Extract main colors
    end
end

% Determine time windows
% Selecting a time window
days = ct{1}; % days of interest 
nights = ct{2}; % nights of interest 
time_window(1) = min([days_crop(days) nights_crop(nights)]);
time_window(2) = max([days_crop(days) nights_crop(nights)]);

% Clean up
clear f p e g t start top_up_bin et ct ac   

%% Parameters Across Time 

% Removing time points with no data
for p = 1:size(parameters,2) % For each parameter
    scrap = []; % Clear scrap
    
    for g = 1:max(group_tags) % For each group
        scrap = [scrap ; find(sum(isnan(parameter_time_sec_smooth{p}(:,group_tags == g)),2) == ...
            size(find(group_tags == g),1))]; % Find time points without data
    end
    
    scrap = unique(scrap); % Clear repeated values
    scrap = sort(scrap); % Sort to ascending order
    
    parameter_time_sec_smooth_y{p}(:,1) = 1:size(parameter_time_sec_smooth{p},1); 
        % Generate a time-line
    parameter_time_sec_smooth_y{p}(scrap,:) = []; % Remove times with no data  
    parameter_time_sec_smooth{p}(scrap,:) = []; % Remove rows with no data 
    
end

clear scrap p g 

%% Distributions

% Calculating the boundaries for each fit 
dist_boundaries = nan(size(wake_cells,2)*2,size(parameters,2),'single');
    % fish*2 (max/min) x parameters 

a = 1; % Start a counter 
for f = 1:size(wake_cells,2) % For each fish  
    dist_boundaries(a:a+1,1:6) = minmax(wake_cells{1,f}(:,3:end)')'; % Wake Bout Parameters 
    dist_boundaries(a:a+1,7:9) = minmax((squeeze(parameter_matrix(f,7:9,:))')')'; % Totals  
    dist_boundaries(a:a+1,10) = minmax(sleep_cells{1,f}(:,3)')'; % Sleep Bout Length  
    dist_boundaries(a:a+1,11:12) = minmax((squeeze(parameter_matrix(f,11:12,:))')')'; % Totals
    a = a + 2; 
end 

% Pre-allocation 
parameter_dists = cell(1,size(parameters,2)); 
for p = find(parameter_smooth == 0) % For most parameters 
    parameter_dists{p} = nan(size(wake_cells,2),size(min(dist_boundaries(:,p)):...
        max(dist_boundaries(:,p)),2),size(parameter_matrix,3),'single'); 
        % {parameters} fish x parameter range x time 
end 

tic
% Distribution Fitting
counter = 0; % Start counter  
for p = find(parameter_smooth == 0) % For most parameters
    
    if ismember(p,1:6) == 1  % For wake parameters
        for f = 1:size(wake_cells,2) % For each fish
            for t = 1:size(parameter_matrix,3) % For each time window                
                
                if isempty(find(parameter_indicies{1,f} == t)) == 0 % If there are bouts
                    pd = fitdist(wake_cells{1,f}(parameter_indicies{1,f}==t,p+2),'kernel','Width',1); % Fit
                    parameter_dists{p}(f,:,t) = pdf(pd,min(dist_boundaries(:,p)):max(dist_boundaries(:,p)));
                else % If there are no bouts
                    parameter_dists{p}(f,:,t) = zeros(1,size(min(dist_boundaries(:,p)):...
                        max(dist_boundaries(:,p)),2)); % Fill with zeros 
                end
            end
        end
              
    elseif p == 10 % For sleep bout length
        for f = 1:size(sleep_cells,2) % For each fish
            for t = 1:size(parameter_matrix,3) % For each time window
                
                if isempty(find(parameter_indicies{2,f} == t)) == 0 % If there are bouts
                    pd = fitdist(sleep_cells{1,f}(parameter_indicies{2,f}==t,3),'kernel','Width',1); % Fit
                    parameter_dists{p}(f,:,t) = pdf(pd,min(dist_boundaries(:,p)):max(dist_boundaries(:,p)));
                    % Note that fitdist ignores NaN values (from water top up)
                else % If there are no bouts 
                    parameter_dists{p}(f,:,t) = zeros(1,size(min(dist_boundaries(:,p)):...
                        max(dist_boundaries(:,p)),2)); % Fill with zeros 
                end
                
            end
        end
        
    end
    
    counter = counter + 1; % Add to counter 
    disp(horzcat('Parameter ',num2str(counter),' of ',...
        num2str(size(find(parameter_smooth == 0),2)),' Fit'));
end
 toc 
 
clear a counter f p pd t time_start time_stop

%% Figure - Parameter Means 

figure; 
for p = 1:size(parameter_comparisons,2) - 2 % For each parameter
    subplot(2,5,p); hold on; clear scrap; counter = 1; 
    set(gca,'FontName','Calibri'); % Set Font  
    title(parameters{p}); % Add title
    
    % Plot parameters
    for e = 1:max(experiment_tags) % For each experiment 
        for g = 1:max(group_tags) % For each group
            clear data; 
            data = parameter_comparisons{p}(:,g,time_window(1):time_window(2)); % extract data 
            data(isnan(data)) = []; % remove nan values 
            data = reshape(data,[group_sizes(g),size(time_window(1):time_window(2),2)]); % reshape to matrix 
            data = data/unit_conversion(1,p); % convert to appropriate units   
            
            legend_lines(g) = errorbar(nanmean(data(experiment_tags(group_tags == g) == e,:)),...
                nanstd(data(experiment_tags(group_tags == g) == e,:))./...
                sqrt(size(find(experiment_tags == e & group_tags == g),1)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);
    
            if e == 1 % For the first experiment 
            legend_cell{g} = horzcat(geno_list.colheaders{g},', n = ',...
                num2str(size(find(group_tags == g),1)));
            % Append the group size to each group name
            end 
            
            % To accurately determine axis scale
            % Add and subtract the std from each mean then determine
            % The highest & lowest value for each group
            scrap(1,counter) = max(nanmean(data(experiment_tags(group_tags == g) == e,:))...
                + nanstd(data(experiment_tags(group_tags == g) == e,:))./...
                sqrt(size(find(experiment_tags == e & group_tags == g),1)));
            scrap(2,counter) = min(nanmean(data(experiment_tags(group_tags == g) == e,:))...
                - nanstd(data(experiment_tags(group_tags == g) == e,:))./...
                sqrt(size(find(experiment_tags == e & group_tags == g),1)));
            
            counter = counter + 1; % Add to counter 
        end
    end 
    
    % Add night patches
    y_lims = [(min(scrap(2,:)) - min(scrap(2,:))*0.05) ...
        (max(scrap(1,:)) + max(scrap(1,:))*0.05)]; % Add a bit of space either side 
        
    a = 1; night_start = first_night; % Start counters  
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
            1 (y_lims(2)-y_lims(1))],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
    
    % Figure Looks
    if p == 5 % For the 4th parameter
        [~,~,~,~] = legend(legend_cell,'Location','best'); % Generate axis
        legend('boxoff'); % Turn legend off  
    end
    axis([0.5 size([days_crop(days) nights_crop(nights)],2)+0.5 ...
        y_lims]); % Set axis 
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format 
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels 
    set(gca, 'XTick', []); % Turn off X-Ticks 
    ylabel(units(p),'Fontsize',12); % Y Labels 

end
  
clear data p scrap g legend_lines y_lims a n r night_start count

%% Figure - Parameters Across Time  

figure; hold on; clear legend_lines
for p = 1:size(parameters,2) % For each parameter
    subplot(3,4,p); hold on; title(parameters{p}); clear scrap; clear y;
    set(gca,'FontName','Calibri'); % Set Font  
    
    y = find(parameter_time_sec_smooth_y{p} >= lb_sec(time_window(1))...
        & parameter_time_sec_smooth_y{p} < lb_sec(time_window(2)+1)); 

    for g = 1:max(group_tags) % For each group
        
        legend_lines(g) = plot(parameter_time_sec_smooth_y{p}(y,1),...
            smooth(nanmean(parameter_time_sec_smooth{p}(y,group_tags == g),2)...
            ,time_bins)/unit_conversion(2,p),'color',cmap(g,:)); % Plot 
                
        legend_cell{g} = horzcat(geno_list.colheaders{g},', n = ',...
            num2str(size(find(group_tags == g),1)));
        % Append the group size to each group name
            
        % Find the top & bottom 
            % Excluding values that smooth outside of the data
        scrap(1,g) = max(legend_lines(g).YData((time_bins):...
            (size(legend_lines(g).YData,2) - time_bins)));
        scrap(2,g) = min(legend_lines(g).YData((time_bins):...
            (size(legend_lines(g).YData,2) - time_bins)));
        
    end
    
    a = 1;
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[lb_sec(nights_crop(nights(n))) ...
            min(scrap(2,:)) - (min(scrap(2,:))*0.05)...
            (lb_sec(nights_crop(nights(n))+1)-1) - lb_sec(nights_crop(nights(n)))...
            max(scrap(1,:)) + (max(scrap(1,:))*0.05)],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1;
    end
    
    if p == 4 % For the 4th parameter - add a legend 
        [~,~,~,~] = legend(legend_cell,'Location','best');
        legend('boxoff'); 
    end
    
    axis([(parameter_time_sec_smooth_y{p}(y(1),1) + time_bins)...
        (parameter_time_sec_smooth_y{p}(y(end),1) - time_bins) ...
        min(scrap(2,:)) - (min(scrap(2,:))*0.05) ...
        max(scrap(1,:)) + (max(scrap(1,:))*0.05)]);
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
    xlabel('Time (Days/Nights)','Fontsize',12);
    set(gca, 'XTick', []);
    ylabel(units_2(p),'Fontsize',12);
end

clear a g legend_lines legend_cell n p r scrap y 

%% Figure - Activity   

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
[~,icons,plots,~] = legend(legend_cols,legend_cell,'Location','best');
legend('boxoff'); 
set(icons(1:max(group_tags)),'Fontsize',32) ; set(plots,'LineWidth',3);
axis([lb_sec(time_window(1)) (lb_sec(time_window(2)+1)-1) 0 y_lims(2)]);
x = lb_sec(time_window(1)):(60*60*12):(lb_sec(time_window(2)+1)-1); 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12})
xlabel('Time (Hours)','Fontsize',32);
ylabel({'Delta Px' ; '(Total/Second)'},'Fontsize',32);

clear a g h icons plots str legend_cell legend_cols legend_lines n r x y_lims 

%% Figure - Parameter Distributions V3
% Uses a log axis for plots
% Plot's Std Rather than SEM

figure;
for p = 1:size(parameters,2) - 2 % For each parameter
    subplot(2,5,p); hold on; col = 1; counter = 1;
    clear data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
    title(parameters{p}); % Add title
    
    if ismember(p,[7:9 11:12]) == 0 % For most parameters
        
        % Take day & night distributions across time windows
        data{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,days_crop(days)),2));
        data{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,nights_crop(nights)),2));
        
        % Expand group tags to account for multiple days/nights
        data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]);
        data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]);
        
        % Expand experiment tags to account for multiple days/nights
        data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
        data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);
        
        % V3
        crop = max(dist_boundaries(:,p)); % crop = all data
        
        for e = 1:max(experiment_tags) % For each experimet
            col = 1;
            for g = 1:max(group_tags) % For each group
                for t = 1:2 % For day/night
                    
                    legend_lines(col) = shadedErrorBar((min(dist_boundaries(:,p)):crop)/unit_conversion(1,p),...
                        nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,:)),...
                        nanstd(data{1,t}(data{2,t} == g & data{3,t} == e,:))...
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
                    y_lims(1,counter) = max(legend_lines(col).mainLine.YData);
                    y_lims(2,counter) = min(legend_lines(col).mainLine.YData);
                    
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
        xlabel(units(p),'Fontsize',12); % X labels
        ylabel('Probability','Fontsize',12); % Y label
        
    else % For the other parameters
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
                    'xValues',g,'distributionColors',...
                    cmap_2(col+1,:)+(1-cmap_2(col+1,:))*(1-(1/e^.5)),'showMM',2);
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
            end
            set(findall(gca,'type','line'),'markersize',15); % change marker sizes

            col = col + 2;
            
        end
        ylabel(units(p),'Fontsize',12); % Y labels
        set(gca,'xtick',1:max(group_tags)); % Set x ticks 
        set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group
    end
    if p == 5 % For the 5th parameter - add a legend
        [~,icons,~,~] = legend(legend_cols,legend_cell,'Location','best');
        legend('boxoff');
        set(icons(1:max(group_tags)*2),'Fontsize',13);
        set(icons((max(group_tags)*2)+1:2:end),'LineWidth',3);
    end
    
end
 
clear col g p spread_cols t xl icons data

%% Figure - Parameter Distributions V3 (With No Error Bars)
% Uses a log axis for plots
% Plot's no Error Bar for ease of Visualization

figure;
for p = 1:size(parameters,2) - 2 % For each parameter
    subplot(2,5,p); hold on; col = 1; counter = 1;
    clear data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
    title(parameters{p}); % Add title
    
    if ismember(p,[7:9 11:12]) == 0 % For most parameters
        
        % Take day & night distributions across time windows
        data{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,days_crop(days)),2));
        data{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,nights_crop(nights)),2));
        
        % Expand group tags to account for multiple days/nights
        data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]);
        data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]);
        
        % Expand experiment tags to account for multiple days/nights
        data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
        data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);
        
        % V3
        crop = max(dist_boundaries(:,p)); % crop = all data
        
        for e = 1:max(experiment_tags) % For each experimet
            col = 1;
            for g = 1:max(group_tags) % For each group
                for t = 1:2 % For day/night
                    
                    legend_lines(col) = plot((min(dist_boundaries(:,p)):crop)/unit_conversion(1,p),...
                        nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,:)),...
                        'color',cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),...
                        'linewidth',3);
                    % Plot - note that this is now in appropriate units
                    % (eg.seconds)
                    
                    legend_cols(col,:) = legend_lines(col).Color; % Save color
                    
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
                    y_lims(1,counter) = max(legend_lines(col).YData);
                    y_lims(2,counter) = min(legend_lines(col).YData);
                    
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
        xlabel(units(p),'Fontsize',12); % X labels
        ylabel('Probability','Fontsize',12); % Y label
        
    else % For the other parameters
        col = 1;
        for g = 1:max(group_tags) % for each group
            clear data;
            % Day
            data{1} = parameter_comparisons{p}(:,g,days_crop(days)); % extract data 
            data{1}(isnan(data{1})) = []; % remove nan values
            data{1} = nanmean(reshape(data{1},[group_sizes(g),size(days_crop(days),2)]),2)/...
                unit_conversion(1,p); % take mean and convert units 
            
            % Night
            data{2} = parameter_comparisons{p}(:,g,nights_crop(nights)); % extract data 
            data{2}(isnan(data{2})) = []; % remove nan values 
            data{2} = nanmean(reshape(data{2},[group_sizes(g),size(nights_crop(nights),2)]),2)/...
                unit_conversion(1,p); % take mean and convert units 
            
            % Plot
            for e = 1:max(experiment_tags)
                spread_cols = plotSpread(data{1}(experiment_tags(group_tags == g) == e),...
                    'xValues',g,'distributionColors',...
                    cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),'showMM',2);
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
                spread_cols = plotSpread(data{2}(experiment_tags(group_tags == g) == e),...
                    'xValues',g,'distributionColors',...
                    cmap_2(col+1,:)+(1-cmap_2(col+1,:))*(1-(1/e^.5)),'showMM',2);
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
            end
            set(findall(gca,'type','line'),'markersize',15); % change marker sizes

            col = col + 2;
            
        end
        ylabel(units(p),'Fontsize',12); % Y labels
        set(gca,'xtick',1:max(group_tags)); % Set x ticks
        set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group
    end
    
    if p == 5 % For the 5th parameter - add a legend
        [~,icons,~,~] = legend(legend_lines,legend_cell,'Location','best');
        legend('boxoff');
        set(icons(1:max(group_tags)*2),'Fontsize',13);
        set(icons((max(group_tags)*2)+1:2:end),'LineWidth',3);
    end
    
end

clear col g p spread_cols t xl icons data counter crop legend_cell legend_cols...
    legend_lines y_lims 

%% Stats - Two Way ANOVA 
    % Note that when fish don't have bouts in a time window, they retain
    % their default NaN value in parameter comparisons, then get filtered 
    % and cause indexing problems 
    % As these cases essentially indicate dead fish I've left the code as
    % it is for now and would suggest excluding these fish in the genotype
    % list
    
% Define groups 
anova_group = []; % group  
anova_experiment = []; % experiments 
for g = 1:max(group_tags) % For each group 
    anova_group = [anova_group  ones(max(group_sizes),1)*g]; % Allocate tag 
    anova_experiment = [anova_experiment ; experiment_tags(group_tags == g)]; % Experiment tags
end 
anova_group = repmat(anova_group(:)',[1,size([days_crop(days) nights_crop(nights)],2)]); % Expand for all days/nights  
anova_experiment = repmat(anova_experiment',[1,size([days_crop(days) nights_crop(nights)],2)]); % Experiment 

anova_time = []; % time 
for t = time_window(1):time_window(2) % For each time window 
    anova_time = [anova_time ones(max(group_sizes),max(group_tags))*mod(t,2)]; 
    % Allocate alternating zeros and ones to each time window 
end 
anova_time = anova_time(:)'; % Vectorise 

if size(days_crop(days),2) == size(nights_crop(nights),2) &&... % If there are an equal number of windows
        size(days_crop(days),2) > 1 % & there is more than one day
    anova_development = []; % development
    anova_development = zeros(1,size(anova_group,2)); % Pre-allocate
    d = 1:size(anova_development,2)/(size(time_window(1):time_window(2),2)/2):...
        size(anova_development,2); % divide into "24h" windows
    for t = 1:size(d,2)-1
        anova_development(d(t):d(t+1)-1) = t;
    end
    
else
    anova_development = ones(1,size(anova_group,2)); % allocate
end

% Calculations
for p = 1:size(parameters,2)-2 % For each parameter
    clear scrap;
    scrap = parameter_comparisons{p}(:,:,time_window(1):time_window(2));
    scrap = scrap(:)'; % Vectorise  
    
    if p == 1 % For the first parameter remove NaN values 
        anova_group(isnan(scrap)) = []; 
        anova_time(isnan(scrap)) = []; 
        if exist('anova_development','var') == 1 % Check if variable exists 
            anova_development(isnan(scrap)) = []; 
        end 
    end 
    
    scrap(isnan(scrap)) = []; % remove NaN values 

    [twa.p(:,p),~,twa.stats{p}] = anovan(scrap,...
        {anova_group,anova_time,anova_development,anova_experiment},...
        'display','off','model','full');
            
end

clear anova_development anova_group anova_time anova_experiment p scrap t

%% Stats - KW 

for p = 1:size(parameters,2)-2 % For each parameter
    clear data data_2; c = 1; d = 1; % Start counters 
    if ismember(p,[7:9 11:12]) == 0 % For most parameters 
        data{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,days_crop(days)),2)); % Day data
        data{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,nights_crop(nights)),2)); % Night data
        
        data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]); % Day tags
        data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]); % Night tags
        
        % Expand experiment tags to account for multiple days/nights
        data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
        data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);
        
        % Compare Day vs Night for each group
        for e = 1:max(experiment_tags) % For each experiment 
            for g = 1:max(group_tags) % For each group
                [kw.time_h(d,p),kw.time_p(d,p)] = kstest2(nanmean(data{1,1}(data{2,1} == g & data{3,1} == e,:)),...
                    nanmean(data{1,2}(data{2,2} == g & data{3,2} == e,:))); % Compare day vs night distributions
                d = d + 1; 
                for n = g+1:max(group_tags) % for each group comparison
                    for t = 1:2 % For each time window
                        [kw.group_h{p}(c,t),kw.group_p{p}(c,t)] = kstest2(nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,:)),...
                            nanmean(data{1,t}(data{2,t} == n & data{3,t} == e,:))); % Make group-wise comparisons
                    end
                    c = c + 1;
                end
            end
        end 
        
    else
        data{1,1} = squeeze(nanmean(permute(parameter_comparisons{p}(:,:,days_crop(days)),[1 3 2]),2)); % Day data
        data{1,2} = squeeze(nanmean(permute(parameter_comparisons{p}(:,:,nights_crop(nights)),[1 3 2]),2)); % Night data
             
        % Split data by group 
        for g = 1:max(group_tags) % for each group 
            data_2{1,g} = data{1,1}(:,g); data_2{1,g}(isnan(data_2{1,g})) = []; % days  
            data_2{2,g} = data{1,2}(:,g); data_2{2,g}(isnan(data_2{2,g})) = []; % nights 
        end 
        
        for e = 1:max(experiment_tags) % For each experiment
            for g = 1:max(group_tags) % For each group
                [kw.time_h(d,p),kw.time_p(d,p)] = kstest2(data_2{1,g}(experiment_tags(group_tags == g) == e,1),...
                    data_2{2,g}(experiment_tags(group_tags == g) == e,1)); % Compare day vs night
                d = d + 1; % add to counter
                for n = g+1:max(group_tags) % for each group comparison
                    for t = 1:2 % For each time window
                        [kw.group_h{p}(c,t),kw.group_p{p}(c,t)] = kstest2(data_2{t,g}(experiment_tags(group_tags == g) == e,1),...
                            data_2{t,n}(experiment_tags(group_tags == n) == e,1)); % Make each comparison
                    end
                    c = c + 1; % add to counter
                end
            end
        end
    end
end

clear p data data_2 c d e g n t 