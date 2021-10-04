%% need-betaine: analysis 1

% goal: intial quantification of cell size distributions from Deanna's 2021-06-25
%       data. are size distributions of mutant vs. wild-type cells
%       different? in normal vs. high salt?

% strategy:
%
% Part ONE: measurements from raw images
%
%   0. initialize experiment data
%   0. define image name of each channel
%   0. for each conditions, build directory (of the 10 images) and loop through each image
%   1. for each image, make a mask to isolate pixels representing cells
%   2. determine which pixels correspond to each cell or "particle"
%   3. measure and record parameters (e.g. size) for each particle
%   4. repeat for all images in all condition
%   5. output from all of this is a data structure "dm" with all particle data,
%      where columns = conditions, rows = an individual image
%
% Part TWO: trim measured data and create final data structure
%
%   1. load "dm" (created in Part One) and define experiment data
%   2. trim measured particules by size
%   3. save new data matrix that groups all data from the same condition
%      into a single data cell
%
% Part THREE:
%
%   6. boxplots to compare width and length across conditions


% ok, let's go!

% last updated: jen, 2021 Sept 29
% commit: 2021-06-25 analysis, first bar plots of width and length


%% Part ONE: measurements from raw images 

clc
clear


% 0. initialize experiment to analyze
experiment = '06-25-21';


% 0. path to raw image files
data_folder = strcat('/Users/jen/Documents/TropiniLab/Data/Deanna/',experiment); 
cd(data_folder)


% 1. initialize experiment meta data
conditions = {'Wt_BHIS','Wt_NaCl1050','Bt1749_BHIS','Bt1749_NaCl1050','Bt1751_BHIS','Bt1751_NaCl1050'};
experiment = '06-25-21';


% 2. initialize image meta data
%px_size = 11/magn; % 11 um pixels with experiment specific magnification
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix



% 3. loop through each condition in experiment and create/add to data matrix, "dm.mat"
%    - each row in "dm.mat" is the meta data or data measured from a given particle.
%    - columns in "dm.mat"
%       1. condition
%       2. image # (as read by MATLAB, not necessarily the suffix in file name)
%       3. width (um)
%       4. length (um)
%       5. area (um^2)
%       6. eccentricity
%       7. orientation
%       8. centroid_x
%       9. centroid_y

for cond = 1:length(conditions)  % length(testset) = # of conditions in experiment
    
    % i. make a directory of images with current condition
    cd(data_folder)
    current_condition = conditions{cond};
    cDirectory = dir(strcat(current_condition,'_*'));
    disp(strcat('Image names of current condition:',current_condition))
    names = {cDirectory.name} % display all names for images of current strain
    
    
    % ii. loop through all images in current condition
    for img = 1:length(names)
        
        % a. determine current image and move to raw data
        cd(data_folder)
        current_img = names{img};
        cd(strcat(current_img,'/Default'))
        
        % b. read in phase image
        img_phase = imread(name_phase); % read phase image
        
        % c. proceed with image processing identical to "whos_a_cell_dp.m"
        phase_smoothed = imgaussfilt(img_phase,0.8); % smooths small pixel noise in background
        bw = edge(phase_smoothed,'sobel'); % detect edges by contrast
        se = strel('disk',1); % structuring element = disk of radius x pixels
        bw_dil = imdilate(bw, se); % dilate
        bw_fill = imfill(bw_dil, 'holes'); % fill
        bw_final = imerode(imerode(bw_fill,se),se); % erode twice

        % d. measure particle parameters
        cc = bwconncomp(bw_final);
        stats = regionprops(cc,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        clear bw bw_dil bw_fill se phase_smoothed  
        
        % e. to data structure, add condition and image data for each particle
        for particle = 1:cc.NumObjects
            stats(particle).condition = current_condition; % condition
            stats(particle).image = img; % image from which that particle came
        end
        
        % f. store particle measurements into cell per image
        dm{img,cond} = stats;
        clear stats
        
    end  
end

cd('/Users/jen/Documents/TropiniLab/Data/Deanna/06-25-21')
save(strcat('dm-raw-',experiment,'.mat'),'dm')


%% Part TWO: trim measured data and create final data structure


clear
clc

% 0. initialize experiment to analyze
experiment = '06-25-21';
magn = 100; % 100x magnification
px_size = 11/magn; % 11 um pixels with
conditions = {'Wt_BHIS','Wt_NaCl1050','Bt1749_BHIS','Bt1749_NaCl1050','Bt1751_BHIS','Bt1751_NaCl1050'};

% 0. initialize segmentation parameters
minWidth = 0.7; % min width for what is a single cell (see whos_a_cell_dp.m)
maxWidth = 3;   % generous max width

% 0. initialize experiment data
cd('/Users/jen/Documents/TropiniLab/Data/Deanna/06-25-21')
load(strcat('dm-raw-',experiment,'.mat'))


 
% 1. concatenate data from same sample
for col = 1:length(conditions)
    
    current_condition = dm(:,col);
    num_imgs = sum(~cellfun(@isempty,current_condition));
    
    combined_cond_data = [];
    for im = 1:num_imgs
        img_data = current_condition{im,1};
        combined_cond_data = [combined_cond_data; img_data];
    end
    clear im img_data
    combined_particles{1,col} = combined_cond_data;
    
end
clear col combined_cond_data num_imgs current_condition img_data



% 2. convert measurements to microns based on imaging specifications     
for jj = 1:length(conditions)
    
    cond_particles = combined_particles{1,jj};
    if isempty(cond_particles) == 1
        continue
    end
    
    
    % 2a. convert x,y coordinate of particle centroid
    x_position = [];
    y_position = [];
    for im = 1:length(cond_particles)
        
        centroid = cond_particles(im).Centroid.*px_size;
        x_position = [x_position; centroid(1)];
        y_position = [y_position; centroid(2)];
        
    end
    parameter_unit.X = x_position;
    parameter_unit.Y = y_position;
    clear particle x_position y_position centroid ii
    
    
    % 2b. convert area
    area = extractfield(cond_particles,'Area')';
    parameter_unit.A = area.*px_size^2;
    clear area
    
    
    % 2c. major & minor axis
    majorAxis = extractfield(cond_particles,'MajorAxisLength')';
    minorAxis = extractfield(cond_particles,'MinorAxisLength')';
    parameter_unit.MajAx = majorAxis.*px_size;
    parameter_unit.MinAx = minorAxis.*px_size;
    clear majorAxis minorAxis
    
    
    % 2d. eccentricity and angle
    ecc = extractfield(cond_particles,'Eccentricity')';
    angle = extractfield(cond_particles,'Orientation')';
    parameter_unit.Ecc = ecc;
    parameter_unit.Angle = angle;
    clear ecc angle
    
    
    % 2e. condition and image
    parameter_unit.Condition = extractfield(cond_particles,'condition')';
    parameter_unit.Image = extractfield(cond_particles,'image')';
    
    
    % 3. trim particles by width
    %    values are set as recorded in whos_a_cell_dp.m
    TrimField = 'MinAx';  % choose relevant characteristic to restrict, run several times to apply for several fields
    p_trim = ParticleTrim_glycogen(parameter_unit,TrimField,minWidth,maxWidth);
    

    % 4. store final data 
    trimmed_data{1,jj} = p_trim;
    
end
clear cond_particles p_trim jj TrimField
    

%% Part THREE: visualize box plots of size data


% 1. gather size data from each condition
widths = cell(1,length(conditions));
lengths = cell(1,length(conditions));
nums = zeros(1,length(conditions));
conds = cell(1,length(conditions));

for kk = 1:length(conditions)

    cond_data = trimmed_data{1,kk};
    widths{kk} = cond_data.MinAx;
    lengths{kk} = cond_data.MajAx;
    %nums(kk) = length(cond_data.MinAx);
    %conds{kk} = ones(length(cond_data.MinAx),1)*kk;
   
end
clear kk cond_data


% 2. format x tick labels
labels = {'Wt_BHIS','Wt_NaCl1050','Bt1749_BHIS','Bt1749_NaCl1050','Bt1751_BHIS','Bt1751_NaCl1050'};


% 3. box plot of cell width across conditions
figure(1)
boxplotGroup(widths,'PrimaryLabels',labels)
xtickangle(90)
ylim([0.2 3])
ylabel('Width (um)')

figure(2)
boxplotGroup(lengths,'PrimaryLabels',labels)
xtickangle(90)
ylim([0 6])
ylabel('Length (um)')


%% Part FOUR: visually affirm box plots

% overlay segmentation data on phase image to affirm validity of
% conclusions in boxplots

clc

% 1. initialize condition and image to analyze
cond2check = 1;
pic = 1;

% 2. initialize image meta data
prefix = 'img_';
suffix = '_position000_time000000000_z000.tif';
name_phase = strcat(prefix,'channel000',suffix);
clear prefix suffix

% 3. isolate image measured data
current = trimmed_data{cond2check};
c_lengths = current.MajAx;
c_pics = current.Image;

i_data = current


%         majorAxes = extractfield(stats,'MajorAxisLength')'.*px_size;
%         minorAxes = extractfield(stats,'MinorAxisLength')'.*px_size;
%         angles = extractfield(stats,'Orientation')';
%         centroid = extractfield(stats,'Centroid');
%         centroid_x = centroid(:,1);
%         centroid_y = centroid(2);


