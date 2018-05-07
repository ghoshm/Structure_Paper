% Vp_Extract 

% marcus.ghosh.11@ucl.ac.uk 

%% Info 

% Extracts data from Viewpoint files with data for every frame

% Input
    % Excel Sheets - output from re-running a viewpoint file  
        % Each with:
        % Rows - every frame for every animal
        % Columns including:
            % Type - data (101) or errors (~=101)
            % Location - which animal each data point came from
            % Data - pixels changed that frame
    % Perl Batch data file (txt) - ouput from running Pearl_Batch_192.m  
    % Genotype list (txt): 
        % Columns - group names  
        % Rows - ROI ID's assigning each fish to a group 
        
% Output 
    % Data structures that can be fed into Vp_Analyse.m  
    % Note that sleep_cells will contain NaN values if the fish water was topped up 
    
%% Assumptions 

% Data comes from 192 regions of interest (ROIS) 
% Experiments start during the day 
% Days are 14hours long, nights are 10hours long 

%% Dependencies 

% Files  
    % dir2
        % https://github.com/ghoshm/Utilities/blob/master/dir2.m 

    % lbmap
        %https://uk.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps

    % Nat Sort Files 
        %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 

    % ProgressBar
        %http://uk.mathworks.com/matlabcentral/fileexchange/6922-progressbar
    
% MATLAB Products 
    % Statistics and Machine Learning Toolbox
    % Curve Fitting Toolbox
    % Parallel Computing Toolbox
    % MATLAB Distributed Computing Server

%% Options
    % For *'s see "Notes on Options" below 
    
% General  
lines_per_sheet = 50000; % Specify the number of data points per Excel sheet 
box = 1; % set which of the two boxes you want to use (*) 
threshold = 200; % Maximum delta px value (**)  
top_up = []; % alter to light boundaries where you topped up fish water. 
    % E.g. Day 2 and 3 (top_up = [2 4]). 
    % E.g. Not topped up (top_up = []). 
time_bins = 60*15; % Choose plot smoothing (in minutes) 
days = [1 2 3 4]; % number of days  
nights = [1 2 3]; % number of nights 

% Figure Settings  
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

% Colors 
col = 'RedBlue'; % (***) 
night_color = [0.9608 0.9608 0.9608]; % color of night shading   

%% Notes on Options  
%* at this point its easier to run each seperatly then merge
    % the data later using Vp_Analyse  
    
% ** A hard threshold is a simple solution, though removing whole bouts
    % including these values is best as it avoids leaving "cut bouts" 
    % From my 3 WT experiments (with lids so no hands), the maximum value I
    % Observe is 165px.
    % From a PTZ Dose response (you may expect higher values) the data
    % doesn't exceed this
    % From looking @ 6 Hcrt experiments and the PTZ experiment - data(data > 0)
    % prctile and std, there is no obvious flexible cut off that could be used. 
    
% *** Choice of 
    % 'Blue'       Single-hue progression to purlish-blue (default)
    % 'BlueGray'   Diverging progression from blue to gray
    % 'BrownBlue'  Orange-white-purple diverging scheme
    % 'RedBlue'    Modified spectral scheme
    
%% Selecting Files  

% Select a folder of Excel Sheets 
folder_path = uigetdir([],'Select a folder of Excel sheets'); % Choose your experiment
folder_open = dir2(folder_path); % Open this folder
disp(horzcat('Running File ',folder_path)); % Report file choice  

% Select a geno_list
[filename, pathname] = uigetfile('*.txt', 'Select a Genotype List'); % Select a geno file
if isequal(filename,0) % If no file is selected
    error('No File Selected') % Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) % Show selected filename
end
geno_list = importdata(strcat(pathname,filename),'\t',2); % Load genolist

clear filename pathname 

% Select a Perl Batch Data File 
[filename, pathname] = uigetfile('*.txt',...
    'Select a Perl Batch Output (DATA) File'); % Select a geno file
if isequal(filename,0) % If no file is selected
    error('No Data File Selected') % Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) % Show selected filename
end

%% Load Data from Excel Sheets 

tic
% Pre-allocation
time = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % time {1}
%pause(30); % Wait for memory 
data_type = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Data type {2} 
%pause(30); % Wait for memory 
fish_id = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Fish id {3}
%pause(30); % Wait for memory 
delta_px = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Delta px {4}
%pause(30); % Wait for memory 
sheet_names = cell(size(folder_open,1),1); % Excel sheet names  

% Ordering by File Name 
for f = 1:size(folder_open,1) % For each excel file 
    sheet_names{f} = folder_open(f).name; % Take it's name 
end % Note that these will be in "computer" order 
    % Ie. 1-10, 100, 1000 etc 

[~,O] = natsortfiles(sheet_names); % Use natsortfiles to sort by file name
    clear sheet_names; % Clear Sheet names 
    
a = 1; % Start a counter 
progress = 0; % Start a timer
data_type_errors = 0; % Start a counter  
order_errors = 0; % Start a counter 
order_errors_size = []; % pre-allocate an empty vector
progressbar('Files') %Initialise progress bars 
for f = O' % For each Excel file
    
    % Load data 
        % Raw data structure 
            % 1 - Time 
            % 2 - Data type 
            % 3 - Fish ID 
            % 4 - Delta px
    fid = fopen(strcat(folder_path,'\',folder_open(f).name)); % Open it
    if f == O(1) % For the first file Skip the first blank lines
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',192 + 1);
    else % For the rest of the files 
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',1);
    end
    
    % Error Handling 
    
    % Dealing with data_type errors (e)
        % Note that this includes cropping the end of the experiment off 
    found = find(raw_text{2} ~= 101); % Find errors
    if isempty(found) ~= 1 % If there are errors 
        for e = 1:size(raw_text,2) % For each data type being imported
            raw_text{e}(found) = []; % Remove errors
        end
        data_type_errors = data_type_errors + 1; % Add to error counter  
    end
    clear found;
    
    % Dealing with ordering errors 
    fish_id_order = str2num(char(raw_text{3})); % Convert str to num
    fish_id_order_check = diff(fish_id_order); % Diff these values
    found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
        isnan(fish_id_order_check) == 0,1,'first'); % Find the first frame drops 
    
    while isempty(found) == 0 % While there are dropped frames 

        % Use a sliding window to find where the pattern normalises
        % Cut backwards
        found_clip_b = [191 1];
        if found - found_clip_b(1) < 1
            found_clip_b = [found + 1 found + 1]; 
        else
            while sum(fish_id_order_check(found - found_clip_b(1):found - found_clip_b(2))) ~= 191
                found_clip_b = found_clip_b + 1; % Slide window
                
                % Catch Running past the start of the file exception
                if found - found_clip_b(1) < 1
                    found_clip_b = [found_clip_b(1) + 1 found_clip_b(1) + 1];
                    break
                end
            end
        end
        
        % Cut forwards
        found_clip_f = [1 191];
        if found + found_clip_f(2) > length(fish_id_order_check)
            found_clip_f = [length(fish_id_order_check) - found + 2 ... 
                length(fish_id_order_check) - found + 2]; 
        else
            while sum(fish_id_order_check(found + found_clip_f(1):found + found_clip_f(2))) ~= 191
                found_clip_f = found_clip_f + 1; % Slide window
                
                % Catch Running past the end of the file exception
                if found + found_clip_f(2) > length(fish_id_order_check)
                    found_clip_f = [found_clip_f(2)+1 found_clip_f(2)+1];
                    break
                end
            end
        end
        
        % Now set values between these sections to NaN
        fish_id_order(found - found_clip_b(2)+2:found + found_clip_f(1)-1) = NaN;
        order_errors_size(order_errors+1,1) = ...
            size(found - found_clip_b(2)+2:found + found_clip_f(1)-1,2); 
            % Store the size of the removed data 
        clear found_clip_b found_clip_f
        
        fish_id_order_check = diff(fish_id_order); % Diff the new fish order
        found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
            isnan(fish_id_order_check) == 0,1,'first'); % Find other frame drops
        order_errors = order_errors + 1; % Add to the order errors counter 
    end
        
    % Data Storage 
    if a == 1 % For the first file 
        time(a:size(raw_text{1},1),1) = raw_text{1}; % Time 
        data_type(a:size(raw_text{1},1),1) = raw_text{2}; % Data Type 
        fish_id(a:size(raw_text{1},1),1) = fish_id_order; % Fish Id order 
        delta_px(a:size(raw_text{1},1),1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter 
    else % For all other files 
        time(a:a+size(raw_text{1},1)-1,1) = raw_text{1}; % Time 
        data_type(a:a+size(raw_text{1},1)-1,1) = raw_text{2}; % Data Type 
        fish_id(a:a+size(raw_text{1},1)-1,1) = fish_id_order; % Fish Id order 
        delta_px(a:a+size(raw_text{1},1)-1,1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter
    end
        
    % Clean up 
    fclose(fid); clear fish_id_order fish_id_order_check found raw_text; % Clear variables
    
    progress = progress + 1; % Add to timer 
    progressbar(progress/size(folder_open,1)); %Update the Fish progressbar
    
end
disp('Time taken to load data from Excel sheets'); % report progress
toc 

% Report back on errors 
disp(horzcat(num2str(data_type_errors),' Files with Data Type Errors')); 
disp(horzcat(num2str(order_errors),' Order Errors')); 
for e = 1:size(order_errors_size,1) % For each order error 
    disp(horzcat('Order Error ',num2str(e),' Size = ',num2str(order_errors_size(e)))); 
end 
disp(horzcat('Ran File ',folder_path)); % Report file choice  

clear a e ans data_type data_type_errors f fid folder_open folder_path ...
    lines_per_sheet O order_errors order_errors_size progress

%% Checking for Persistant Ordering Errors
% At this point only one type of ordering error will remain
% This will be the rare case where the ordering error bridges two Excel sheets
% E.g. 171030_16

% Check for remaining Errors
scrap = diff(isnan(fish_id));
% This will leave either
% +1 (A number to a NaN)
% -1 (A NaN to a number)
% 0 (number to number)

% Check for Backwards errors
check(:,1) = find(scrap == 1); % Locs
check(:,2) = fish_id(scrap == 1); % Values (should all be 192)

if isempty(find(check(:,2) ~= 192)) == 0 % If there are errors
    err = find(check(:,2) ~= 192); % define them
    fish_id_diff = diff(fish_id); % diff fish_id
    
    for e = 1:size(err,1) % for each error
        
        % Use a sliding window to find where the pattern normalises
        % Cut backwards
        found = check(err(e,1),1); % set the error location
        found_clip_b = [191 1];
        if found - found_clip_b(1) < 1
            found_clip_b = [found + 1 found + 1];
        else
            while sum(fish_id_diff(found - found_clip_b(1):found - found_clip_b(2))) ~= 191
                found_clip_b = found_clip_b + 1; % Slide window
                
                % Catch Running past the start of the file exception
                if found - found_clip_b(1) < 1
                    found_clip_b = [found_clip_b(1) + 1 found_clip_b(1) + 1];
                    break
                end
            end
        end
        
        fish_id(found - found_clip_b(2)+2:found) = NaN;
    end
    
end

clear check err e found found_clip_b fish_id_diff

% Check for forwards errors
check(:,1) = find(scrap == -1) + 1;
check(:,2) = fish_id(find(scrap == -1)+1);

if isempty(find(check(:,2) ~= 1)) == 0 % If there are errors
    err = find(check(:,2) ~= 1); % define them
    fish_id_diff = diff(fish_id); % diff fish_id
    
    for e = 1:size(err,1) % for each error
        
        % Cut forwards
        found = check(err(e,1),1); % set the error location
        found_clip_f = [1 191];
        if found + found_clip_f(2) > length(fish_id_diff)
            found_clip_f = [length(fish_id_diff) - found + 2 ...
                length(fish_id_diff) - found + 2];
        else
            while sum(fish_id_diff(found + found_clip_f(1):found + found_clip_f(2))) ~= 191
                found_clip_f = found_clip_f + 1; % Slide window
                
                % Catch Running past the end of the file exception
                if found + found_clip_f(2) > length(fish_id_diff)
                    found_clip_f = [found_clip_f(2)+1 found_clip_f(2)+1];
                    break
                end
            end
        end
        
        fish_id(found:found + found_clip_f(1)-1) = NaN;
    end
    
end

clear check err fish_id_diff e found found_clip_f fish_id_diff scrap

disp('Found & Corrected Persistant Ordering Errors'); % Report

%% Reshape The Data 

% Check that the tracking start's with fish 1 
if fish_id(1) ~= 1 % if not 
   crop = find(fish_id == 1,1,'first') - 1; % find just before fish 1   
   % Crop all of the data 
       delta_px(1:crop) = []; 
       fish_id(1:crop) = []; 
       time(1:crop) = []; 
   clear crop; 
end 

% Check the number of frames per fish is correct 
frames_per_fish = zeros(1,max(fish_id)); 
for f = 1:max(fish_id) % For each fish 
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    frames_per_fish(f) = size(found,1); % Store number of frames 
    disp(horzcat('Found frames for fish Number ',num2str(f),' of ',...
        num2str(max(fish_id)))); % Report on progress 
end 

if min(frames_per_fish) ~= max(frames_per_fish) % If the number of frames is not equal 
   error('Data formatted Incorrectly'); % Call an error  
end 

% Delta px sq 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        delta_px_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    delta_px_sq(:,f) = delta_px(found); % Take delta_px values
    disp(horzcat('Collected frames for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end 
clear f found delta_px

% Time sq 
    % Note this is separated from delta px sq for ease of memory handling 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        time_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    time_sq(:,f) = time(found); % Take time values 
    disp(horzcat('Collected time for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end
clear f found time fish_id frames_per_fish

%% Add in time 

fid = fopen(strcat(pathname,filename)); % Open time data 
formatSpec = '%*884s%f%[^\n\r]';
dataArray = textscan(fid, formatSpec, 3-3+1, 'Delimiter',...
    '', 'WhiteSpace', '', 'HeaderLines', 2, 'ReturnOnError', false, 'EndOfLine', '\r\n');
start_time = dataArray{1}; % Extract start time (Matlab generated code) 
fclose(fid); % Close it 

clear ans dataArray fid filename formatSpec pathname

% Light Boundaries Calculation
time_sq_max = ((max(time_sq'))/(60*60))'; % Max time @ each frame (hours from start) 

a = 1; % start a counter
time_counter = 0; % Start a counter 
boundary = 14 - start_time; % Assumes the experiment starts during the day
while time_counter < time_sq_max(end) - 10 % (10 allows for at least a night)  
    lb(a,1) = knnsearch(time_sq_max,boundary); % Find the best match 
    if mod(a,2) == 1 % If odd
        boundary = boundary + 10; % Add night hours
    else
        boundary = boundary + 14; % Add day hours
    end
    time_counter = time_sq_max(lb(a,1)); % Set time counter 
    disp(horzcat('Found Light Boundary = ',num2str(a))); % Report progress
    a = a + 1; % Add to counter
end

% Day vs Night
dn = ones(size(time_sq,1),1,'single'); % Pre-allocate 
for t = 1:2:size(lb,1) % For each night boundary
    dn(lb(t):lb(t+1)-1) = 0;   
end 

clear a boundary start_time time_counter time time_sq t 

%% Organise the data 

if box == 1
    delta_px_sq(:,97:end) = []; % Remove the unused box
else
    delta_px_sq(:,1:96) = []; % Remove the unused box  
end

delta_px_sq = delta_px_sq - 1; % Set minimum value to zero 

%% Remove "Noise" - set values to zero
% 1. Abnormally high viewpoint values
% 2. Topping up Fish Water

% 1. Remove High Viewpoint values
% Note that this is adapted from the parameter extraction code (below)
wake_cells = cell(1,size(delta_px_sq,2)); % Wake Cells (bout parameters)

% Finding transitions
delta_px_sq_scrap = delta_px_sq;
delta_px_sq_scrap(delta_px_sq_scrap > 0) = 1; % Find active frames
delta_px_sq_scrap = diff(delta_px_sq_scrap); % Diff to find transitions
% 1 = inactive to active
% -1 = active to inactive

for f = 1:size(delta_px_sq,2) % For each fish
    % Note this this runs apporximately twice as fast as just using a
    % For loop
    
    % Starts - ensures no bouts are lost at the start
    if  delta_px_sq(1,f) > 0 % If active in first bin
        wake_cells{1,f}(:,1) = [1 ; find(delta_px_sq_scrap(:,f) == 1)+1]; % Find active bout starts
    else % Ie. if inactive in first bin
        wake_cells{1,f}(:,1) = find(delta_px_sq_scrap(:,f) == 1)+1; % Find active bout starts
    end
    
    % Ends - ensures no bouts are lost at the end
    if delta_px_sq(size(delta_px_sq,1),f) > 0 % If active in last bin
        wake_cells{1,f}(:,2) = [find(delta_px_sq_scrap(:,f) == - 1);...
            size(delta_px_sq,1)]; % Find active bout ends
    else
        wake_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == - 1);
    end
    
    % Parameter extraction
    wake_cells{1,f}(:,3) = NaN; % Pre-allocate
    
    % Active bouts
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        wake_cells{1,f}(b,3) = nanmax(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Max
        
        if wake_cells{1,f}(b,3) > threshold % Hard coded threshold
            delta_px_sq(wake_cells{1,f}(b,1):...
                wake_cells{1,f}(b,2),f) = 0; % Bin to zero
        end
        
    end
    
end

clear b delta_px_sq_scrap f wake_cells threshold

% 2. Filter out Hands & Truncated Bouts - V2 
%figure; plot(nanmax(delta_px_sq')); title('Max'); 
%figure; plot(nanmean(delta_px_sq')); title('Mean'); 

if isempty(top_up) == 0 % if fish h20 was topped up
    top_up_bin = nan(size(top_up,2),2,'single'); % pre-allocate (top ups x start/stop)
    
    fps = round(1/((nanmean(diff(time_sq_max)))*(60*60)),1); % Calculate frame rate
    
    for t = 1:size(top_up,2) % For each top up
        [~,top_up_bin(t,1)] = find(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)'))) >= ...
            nanmean(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))) + ...
            10*nanstd(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))),1,'first'); % find start 
        
        [~,top_up_bin(t,2)] = find(abs(diff(max(delta_px_sq(lb(top_up(t)):(lb(top_up(t)) + top_up_bin(t,1) + (fps*60*20)),:)'))) >= ...
            nanmean(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))) + ...
            10*nanstd(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))),1,'last'); % find stop 
        % Note that here I only look within 20mins (this helps avoid
        % reminaing viewpoint glitches "(fps*60*20)") 
        
        %figure; hold on;
        %plot(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)'))));
        %plot([top_up_bin(t,1) top_up_bin(t,1)] - (fps*90),[0 200],'r','linewidth',3);
        %plot([top_up_bin(t,2) top_up_bin(t,2)] + (fps*90),[0 200],'r','linewidth',3);
        
        % To ensure you get all of the noise, cut a bit more either side 
        top_up_bin(t,1) = lb(top_up(t)) + top_up_bin(t,1) - (fps*90); % go 90s further back
        top_up_bin(t,2) = lb(top_up(t)) + top_up_bin(t,2) + (fps*90); % go 90s further forwards
        
        for f = 1:size(delta_px_sq,2) % for each fish
            if delta_px_sq(top_up_bin(t,1),f) == 0 && delta_px_sq(top_up_bin(t,2),f) == 0
                delta_px_sq(top_up_bin(t,1):top_up_bin(t,2),f) = 0; % set these values to zero
            else % if they have bouts overlapping with these cuts
                delta_px_sq(top_up_bin(t,1)- ...
                    (find(flip(delta_px_sq(1:top_up_bin(t,1),f)) == 0,1,'first')-2):...
                    top_up_bin(t,2) + (find(delta_px_sq(top_up_bin(t,2):end,f) == 0,1,'first')-2),f) = 0;
            end
        end
        

        
    end
    
else % if not topped up
    top_up_bin = []; % store a blank variable
end

clear f fps t top_up

%figure; plot(nanmax(delta_px_sq')); title('Max - Post '); 
%figure; plot(nanmean(delta_px_sq')); title('Mean - Post '); 

%% Group the data by condition 

% Generate group tags 
group_tags = nan(size(delta_px_sq,2),1); % Pre-allocate
for g = 1:size(geno_list.data,2) % For each group 
    group_tags(geno_list.data(1:find(isnan(geno_list.data(:,g))==0,1,'last')...
        ,g)) = g; % Assign group membership  
end 
delta_px_sq(:,isnan(group_tags)) = []; % Remove data   
group_tags(isnan(group_tags)) = []; % Remove blank values 

clear g  

%% Extract Parameters from Data (< Re-shaped data)

% Variables 
tic
% Calculate an approximate frame rate 
fps = round(1/((nanmean(diff(time_sq_max)))*(60*60)),1); 
    % Diff time_sq_max to find time between frames 
    % Take a mean, Convert to mins then hours 
    % Divide one by this value and round -> frames per second
    
parameters = {'Active Bout Length','Active Bout Mean',...
    'Active Bout Variance','Active Bout Total',...
    'Active Bout Minimum','Active Bout Maximum','Number of Active Bouts',...
    'Total Time Active','Total Activity','Inactive Bout Length','Number of Inactive Bouts',...
    'Total Time Inactive'}; % Specify parameters  

% Specify how to convert frame values to parameter units (for figures)
    % Note that these are denominators (will later be used for division) 
unit_conversion(1,:) = [fps 1 1 1 1 1 1 (fps*3600) 1 fps 1 (fps*3600)];
unit_conversion(2,:) = [fps 1 1 1 1 1 1 1 1 fps 1 1];  

% Specify units (for figures) 
units = {'Seconds','Delta Px','Delta Px','Delta Px','Delta Px','Delta Px',...
    'No.','Hours','Delta Px','Seconds','No.','Hours'};
units_2 = {'Seconds','Delta Px','Delta Px','Delta Px','Delta Px','Delta Px',...
    'No.',horzcat('Seconds/',num2str(fps),'s'),'Delta Px','Seconds',...
    'No.',horzcat('Seconds/',num2str(fps),'s')};

% Specify Smoothing operation (0 = mean, 1 = total, 2 = max) - for figures 
parameter_smooth(1:size(parameters,2)) = 0; 
parameter_smooth(7:8) = 1; parameter_smooth(9) = 2; parameter_smooth(11:12) = 1; 

% Pre-allocate 
wake_cells = cell(1,size(delta_px_sq,2)); % Wake Cells (bout parameters)  
sleep_cells = cell(1,size(delta_px_sq,2)); % Sleep Cells (bout parameters) 

for p = 1:size(parameters,2) % For each parameter
    parameter_time{p} = nan(size(delta_px_sq),'single'); % Parameters across time 
end 

% Finding transitions
delta_px_sq_scrap = delta_px_sq; 
delta_px_sq_scrap(delta_px_sq_scrap > 0) = 1; % Find active frames 
delta_px_sq_scrap = diff(delta_px_sq_scrap); % Diff to find transitions  
    % 1 = inactive to active 
    % -1 = active to inactive 
 
parfor f = 1:size(delta_px_sq,2) % For each fish 
        % Note this this runs apporximately twice as fast as just using a
        % For loop 
    
    % Starts - ensures no bouts are lost at the start  
    if  delta_px_sq(1,f) > 0 % If active in first bin  
        wake_cells{1,f}(:,1) = [1 ; find(delta_px_sq_scrap(:,f) == 1)+1]; % Find active bout starts
        sleep_cells{1,f}(:,1) = (find(delta_px_sq_scrap(:,f) == -1)+1); % Find sleep bout starts 
    else % Ie. if inactive in first bin 
        wake_cells{1,f}(:,1) = find(delta_px_sq_scrap(:,f) == 1)+1; % Find active bout starts 
        sleep_cells{1,f}(:,1) = [1 ; (find(delta_px_sq_scrap(:,f) == -1)+1)]; % Find sleep bout starts 
    end 
    
    % Ends - ensures no bouts are lost at the end 
    if delta_px_sq(size(delta_px_sq,1),f) > 0 % If active in last bin 
        wake_cells{1,f}(:,2) = [find(delta_px_sq_scrap(:,f) == - 1);...
            size(delta_px_sq,1)]; % Find active bout ends
        sleep_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == 1); % Find sleep bout ends
    else 
        wake_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == - 1); 
        sleep_cells{1,f}(:,2) = [(find(delta_px_sq_scrap(:,f) == 1)) ; size(delta_px_sq,1)]; % Find sleep bout ends
    end
    
    % Parameter extraction 
    wake_cells{1,f}(:,3:8) = NaN; % Pre-allocate 
    wake_cells{1,f}(:,3) = (wake_cells{1,f}(:,2)+1) - wake_cells{1,f}(:,1); % Wake Bout Length 
    sleep_cells{1,f}(:,3) = (sleep_cells{1,f}(:,2)+1) - sleep_cells{1,f}(:,1); % Sleep Bout Length 

    % Removing "Hands" - setting length to NaN
    try % "check" if the top up variable exists  
        for t = 1:size(top_up_bin,1) % for each top up
            temp = sleep_cells{1,f}(sleep_cells{1,f}(:,3) >= ...
                diff(top_up_bin(t,:)),1:2); % filter for bouts long enough
            idx = knnsearch(temp,top_up_bin(t,:)); % find the best match among filtered bouts
            idx = knnsearch(sleep_cells{1,f}(:,1:2), temp(idx,:)); % find this in the full set
            sleep_cells{1,f}(idx,3) = NaN; % set length to NaN
        end
    catch
    end
    
    % Active bouts 
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        wake_cells{1,f}(b,4) = nanmean(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Mean
        wake_cells{1,f}(b,5) = nanvar(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Variance
        wake_cells{1,f}(b,6) = nansum(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Total
        wake_cells{1,f}(b,7) = nanmin(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Min
        wake_cells{1,f}(b,8) = nanmax(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Max
    end
    
end 
 
% Parameter Time 
for f = 1:size(delta_px_sq,2) % For each fish
    
    % Active Bouts 
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        parameter_time{1}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,3); % Fill in bout length
        parameter_time{2}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,4); % Fill in bout Mean
        parameter_time{3}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,5); % Fill in bout Variance 
        parameter_time{4}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,6); % Fill in bout Total  
        parameter_time{5}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,7); % Fill in bout Minimum 
        parameter_time{6}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,8); % Fill in bout Maximum 
        parameter_time{7}(wake_cells{1,f}(b,1),f) = 1; % No. of Active bouts 
        parameter_time{8}(wake_cells{1,f}(b,1):wake_cells{1,f}(b,2),f) = 1; % Total Time Active  
        parameter_time{9}(wake_cells{1,f}(b,1):wake_cells{1,f}(b,2),f) = ...
            nansum(wake_cells{1,f}(1:b,6)); % Total Activity
    end
    
    % Inactive Bouts 
    for b = 1:size(sleep_cells{1,f},1) % For each sleep bout
        if isnan(sleep_cells{1,f}(b,3)) == 0 % check it's not a "hand artefact"
            parameter_time{10}(sleep_cells{1,f}(b,1),f) = ...
                sleep_cells{1,f}(b,3); % Fill in bout length
            parameter_time{11}(sleep_cells{1,f}(b,1),f) = 1; % No. of Inactive Bouts
            parameter_time{12}(sleep_cells{1,f}(b,1):sleep_cells{1,f}(b,2),f) = 1; % Total Time Inactive
        end
    end
    
    disp(horzcat('Calculated parameters across time for fish = ',num2str(f),...
    ' of ',num2str(size(delta_px_sq,2)))); % Report progress 
end

toc
clear b delta_px_sq_scrap f p   

%% Statistics & Plots - Variables

% Determine day/night order
    % Note that dn currently assumes the experiment starts during the day
lb = [1 ; lb]; % Add 1 to lb
if dn(1) == 1 % If the experiment starts in the day 
    lb_days = lb(1:2:size(lb,1)); % Assign day start values (in frames) 
    lb_nights = lb(2:2:size(lb,1)); % Assign night start values (in frames) 
    days_crop = 1:2:size(lb,1); nights_crop = 2:2:size(lb,1); 
            % Assign logical indicies  
else
    lb_days = lb(2:2:size(lb,1)); % Assign day start values (in frames) 
    lb_nights = lb(1:2:size(lb,1)); % Assign night start values (in frames) 
    days_crop = 2:2:size(lb,1); nights_crop = 1:2:size(lb,1); 
            % Assign logical indicies  
end 

% Determine which windows are of interest 
time_window(1) = min([days_crop(days) nights_crop(nights)]);  
time_window(2) = max([days_crop(days) nights_crop(nights)]); 

% Determine first night  
if min(days_crop(days)) < min(nights_crop(nights))
    first_night = 2; 
else 
    first_night = 1; 
end

% Colours 
cmap_2 = lbmap(max(group_tags)*2,col); % Generate a 2x color map
cmap = cmap_2(1:2:size(cmap_2,1),:); % Extract main colors

% Group sizes 
for g = 1:max(group_tags) % For each group 
    group_sizes(g) = size(find(group_tags == g),1); 
end 

clear g

%% Parameters - Generating Averages (in frames) 
parameter_matrix = nan(size(wake_cells,2),size(parameters,2),...
    size(lb,1)); % Fish x parameters x time windows
parameter_indicies = cell(2,size(wake_cells,2)); % wake/sleep x fish
lb = [lb ; size(delta_px_sq,1)]; % Add end to lb

for f = 1:size(wake_cells,2) % For each fish
    for t = 1:size(parameter_matrix,3) % For each time window
        % Wake bouts
        clear time_start time_stop;
        % Find the first bout that starts within the window
        time_start = find(wake_cells{1,f}(:,1) >= lb(t),1,'first');
        % Find the last bout that starts within the window
        if t+1 < size(lb,1) % For most windows
            time_stop = find(wake_cells{1,f}(:,1) < lb(t+1),1,'last');
        else % For the last window
            time_stop = find(wake_cells{1,f}(:,1) <= lb(t+1),1,'last');
        end
        
        % Store logical index
        parameter_indicies{1,f} = [parameter_indicies{1,f} ; ...
            ones(size(time_start:time_stop,2),1)*t];
        
        % Extract bout parameters (1-6)
        parameter_matrix(f,1:(size(wake_cells{1,f},2)-2),t)...
            = nanmean(wake_cells{1,f}(time_start:time_stop,3:end)); % Means
        % Number of bouts (7)
        parameter_matrix(f,7,t) = size(wake_cells{1,f}...
            (time_start:time_stop,3:end),1);
        % Total time active (8) - sum of lengths
        parameter_matrix(f,8,t) = nansum(wake_cells{1,f}...
            (time_start:time_stop,3),1);
        % Total activity (9) - sum of activity
        parameter_matrix(f,9,t) = nansum(wake_cells{1,f}...
            (time_start:time_stop,6),1);
        
        % sleep bouts (10-12)
        clear time_start time_stop;
        % Find the first bout that starts within the window
        time_start = find(sleep_cells{1,f}(:,1) >= lb(t),1,'first');
        % Find the last bout that starts within the window
        if t+1 < size(lb,1) % For most windows
            time_stop = find(sleep_cells{1,f}(:,1) < lb(t+1),1,'last');
        else % For the last window
            time_stop = find(sleep_cells{1,f}(:,1) <= lb(t+1),1,'last');
        end
        
        % Store logical index
        parameter_indicies{2,f} = [parameter_indicies{2,f} ; ...
            ones(size(time_start:time_stop,2),1)*t];
        
        % Sleep Bout Length (10)
        parameter_matrix(f,10,t)...
            = nanmean(sleep_cells{1,f}(time_start:time_stop,3));
        % Number of bouts (11)
            % Subtract number of NaN's (H20 Top up)
        parameter_matrix(f,11,t) = size(sleep_cells{1,f}...
            (time_start:time_stop,3:end),1)... 
            - sum(isnan(sleep_cells{1,f}(time_start:time_stop,3))); 
        % Total time inactive (12) - sum of lengths
        parameter_matrix(f,12,t) = nansum(sleep_cells{1,f}...
            (time_start:time_stop,3),1);
        
    end
end

% Re-group data for ease of comparisons
parameter_comparisons = cell(1,size(parameters,2)); % Pre-allocate
for p = 1:size(parameter_comparisons,2) % For each parameter
    parameter_comparisons{p}(1:max(group_sizes),1:max(group_tags),...
        1:size(parameter_matrix,3)) = NaN; % {parameters} Most fish per group x
    % x each group x time windows
end

for p = 1:size(parameter_comparisons,2) % For each parameter
    for g = 1:max(group_tags) % For each group
        for t = 1:size(parameter_matrix,3) % For each time window
            parameter_comparisons{p}(1:group_sizes(g),g,t) = ...
                parameter_matrix(group_tags == g,p,t);
        end
    end
end

clear f t time_start time_stop p g

%% Smoothing data into seconds

% Pre-allocate
delta_px_sq_sec = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish
delta_px_sq_sec_smooth = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish
dn_sec = nan(size(1:(fps):size(delta_px_sq,1),2),1,'single'); % time x 1
for p = 1:size(parameters,2) % For each parameter
    parameter_time_sec_smooth{p} = nan(size(delta_px_sq_sec_smooth),'single');
    % {parameters} time (seconds) x fish
end

% Smooth each fish's activity & parameters
for f = 1:size(delta_px_sq,2) % For each fish
    
    a = 1; % Start a counter
    for t = 1:(fps):size(delta_px_sq,1) % For each second
        if t + (fps-1) < size(delta_px_sq,1) % Check to prevent running off the end
            delta_px_sq_sec(a,f) = nansum(delta_px_sq(t:t+(fps-1),f)); % Bin activity
            
            for p = 1:size(parameters,2) % For each parameter
                if parameter_smooth(p) == 0 % For most parameters
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nanmean(parameter_time{p}(t:t+(fps-1),f)); % Mean
                elseif parameter_smooth(p) == 1
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nansum(parameter_time{p}(t:t+(fps-1),f)); % Sum
                elseif parameter_smooth(p) == 2
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nanmax(parameter_time{p}(t:t+(fps-1),f)); % Max
                end
            end
            
            if f == 1 % For the first fish
                dn_sec(a,1) = mode(dn(t:t+(fps-1),1));
                % Take the most common light value within each bin
                % Will need to create indicies for each time window here
                % (04.08.17)
            end
        
        a = a + 1; % Add to counter
        
        else 
            delta_px_sq_sec(a,f) = 0; % This prevents smoothing errors
            dn_sec(a,1) = dn_sec(a-1,1); % Assume last value 
        end
        
    end
    
    % Smooth the activity data
    delta_px_sq_sec_smooth(:,f) = smooth(delta_px_sq_sec(:,f),time_bins);
    disp(horzcat('Smoothed fish ',num2str(f),' of ' ,...
        num2str(size(delta_px_sq,2)))); % Report progress
end

% Determine Day/Night Transitions in seconds
% Note that the binning will make this not 100% accurate
lb_sec = [1 ; find(diff(dn_sec) ~= 0) + 1; size(delta_px_sq_sec_smooth,1)];

clear a f p t 