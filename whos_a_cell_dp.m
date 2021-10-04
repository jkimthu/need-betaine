%% visualize ID'ed particles to determine trimming parameters


% goal: measure cell size from a phase image,
%       overlay measurement onto phase image,
%       visualize and repeat to determine desired values
%       (e.g. values that identify single cells)


% strategy:
%
%   For each condition of interest, user is prompted to:
%       1. choose one test image to identify single cells based on width
%             - have image name ready to choose its index from displayed list
%       2. user is prompted to enter guesses of min and mix width to define
%          a single cell
%       3. user then checks particle overlaps onto phase image to determine
%           if input widths are accurate
%       4. if not accurate, user can adjust guesses and go again
%       5. if acceptable, user can then choose to save data
%       6. this repeats for all unique strain samples (fixed cultures)


% last updated: jen, 2021 Aug 16
% commit: created segdata using same parameters across conditions


% ok, let's go!

%% Part 0. intitialize user specific file paths & experiment of interest

clc
clear

% 1. path to data files and saved segmentation parameters
cd('/Users/jen/Documents/TropiniLab/Data/Deanna/')

% 2. load stored data and define experiment of interest
load('segdata_dp.mat')
tempdata = segdata_dp;

%% Part 1. collect and store segmentation parameters for a new experiment

% 1. initialize experiment meta data
magn = 100;
conditions = {'Wt_BHIS','Wt_NaCl1050','Bt1749_BHIS','Bt1749_NaCl1050','Bt1751_BHIS','Bt1751_NaCl1050'};
experiment = '06-25-21';

% 2. initialize image meta data
px_size = 11/magn; % 11 um pixels with experiment specific magnification
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix

% 3. access image data
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Deanna/',experiment); 
cd(data_folder)


% 4b. for each unique strain, input index chosen image in list of image names
for condition = 1:length(conditions)
    
    % 4b,i. make a directory of images with current condition
    current_condition = conditions{condition};
    cDirectory = dir(strcat(current_condition,'_*'));
    disp(strcat('Image names of current condition:',current_condition))
    names = {cDirectory.name} % display all names for images of current strain
    
    % 4b,ii. prompt user for image to analyse
    disp(strcat('Current condition = ',current_condition))
    prompt = 'Enter one image index in "names" as a double: ';
    img2test = input(prompt);
    condition_img = names{img2test}; % chosen image
    testset{condition} = condition_img;
end
clear strain_groups ss prompt
clear current_condition current_sample strain_img



% 5. loop through each condition in experiment to find segmentation parameters
%    and store parameters in segdata_dp.mat
counter = 1;
while counter <= length(testset)  % length(testset) = # of conditions in experiment
    
    cond = conditions{counter};
    cond_img = testset{counter};
    
    % 5a. display phase image & mask
    cd(data_folder)
    cd(strcat(cond_img,'/Default'))
    img_phase = imread(name_phase); % read phase image
    figure(1) % display phase image
    imshow(img_phase, 'DisplayRange',[1000 5000]); %lowering right # increases num sat'd pxls
    title('phase')
    
    phase_smoothed = imgaussfilt(img_phase,0.8); % smooths small pixel noise in background
    figure(2) % display 
    imshow(phase_smoothed, 'DisplayRange',[1000 5000]);
    
    bw = edge(phase_smoothed,'sobel'); % detect edges by contrast
    figure(3) % display 
    imshow(bw)
    
    se = strel('disk',1); % structuring element = disk of radius x pixels
    bw_dil = imdilate(bw, se); % dilate
    figure(4)
    imshow(bw_dil)
    
    bw_fill = imfill(bw_dil, 'holes'); % fill
    figure(5)
    imshow(bw_fill)
    
    bw_final = imerode(imerode(bw_fill,se),se); % erode twice
    figure(6)
    imshow(bw_final)
    title(strcat('mask'))
    
    figure(7) % display phase image
    imshow(img_phase, 'DisplayRange',[1000 5000]); %lowering right # increases num sat'd pxls
    title(strcat('overlay'))
    
    % 5b. measure particle parameters
    cc = bwconncomp(bw_final);
    stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
    majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
    minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
    angles = extractfield(stats,'Orientation')';
    clear bw bw_dil bw_fill se phase_smoothed
    
    % 5c. prompt user to input min and max width thresholds for single cells
    disp(strcat('Current condition = ',cond))
    prompt_min = 'Enter minimum width (um) to define a single cell as a double: ';
    minWidth = input(prompt_min);
    prompt_max = 'Enter maximum width (um) to define a single cell as a double: ';
    maxWidth = input(prompt_max);
    clear prompt_min prompt_max
    
    % 5d. overlay particle outlines with blue indicating particles within width thresholds
    for pp = 1:length(stats)
        
        centroid = stats(pp).Centroid.*px_size;
        [x_rotated, y_rotated] = drawEllipse2(pp, majorAxes, minorAxes, centroid(1), centroid(2), angles, px_size);
        lineVal = 1;
        
        if minorAxes(pp) < minWidth  % noise threshold
            color = rgb('HotPink');
            color_text = rgb('DarkMagenta');
        elseif minorAxes(pp) > maxWidth % clump threshold
            color = rgb('DeepPink');
            color_text = rgb('DarkMagenta');
        else
            color = rgb('DarkTurquoise'); % single cells
            color_text = rgb('MidnightBlue');
        end
        
        figure(7)
        hold on
        plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
        text((centroid(1)+1)/px_size, (centroid(2)+1)/px_size, num2str(pp),'Color',color_text,'FontSize',10)
        
    end
    
    
    % 5d. prompt user to accept or edit thresholds
    %     i) if accept, store thresholds and proceed to next unique strain
    %    ii) if edit, repeat with current strain until satisfied
    disp(strcat('Current strain sample = ',cond))
    prompt_accept = 'Satisfied? Input "Yes" (stores seg data) or "No" (repeats viz) as a string: ';
    isSatified = input(prompt_accept);
    if strcmp(isSatified, "Yes") == 1
        disp('Saving min and max widths to segdata.mat')
        
        % saving data to segdata.mat
        % a. determine which row to add new data 
        counter = counter + 1;
        newrow = length(tempdata) + 1;
        newdata(1).condition = cond;
        newdata(1).tested_img = cond_img;
        newdata(1).test_experiment = experiment;
        newdata(1).minWidth = minWidth;
        newdata(1).maxWidth = maxWidth;
        
        % b. store data into temporary structure
        tempdata{newrow,1} = newdata;
        
    elseif strcmp(isSatified, "No") == 1
        disp('Thanks for being thorough! Please try again, adjusting min and max widths')
    else
        disp('Error: must input either "Yes" or "No", please restart :P')
    end
    close(figure(1),figure(2),figure(3),figure(4),figure(5),figure(6),figure(7))
    
end

tempdata
segdata_dp

%% Part 2. if happy, save data!
cd('/Users/jen/Documents/TropiniLab/Data/Deanna/')
segdata_dp = tempdata;
save('segdata_dp.mat','segdata_dp')

