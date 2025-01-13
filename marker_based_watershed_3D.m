%% clears command window, clears variables, closes all figures

clc
clear all
close all

%% 1) Load .tiff images into Matlab

% pixel information
pixell=0.5979761; % microns per pixel
pixelw=pixell; % microns per pixel
pixelarea=pixell*pixelw; % microns squared

% stepsize
stepsize=2.0; %microns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% read in NUC images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change current directory to NUC folder
cd('./NUC')

% file extension
file_extension='*.png'; % searches for .tiff files
file_info=dir(file_extension);

% algorithm to delete pseudo-files
indx_del = []; 
for k = 1 : length(file_info)
    filename=file_info(k).name;
    if startsWith(filename, '._')
        indx_del(k) = k;
    end
end

if exist('indx_del','var') == 1
    indx_del = nonzeros(indx_del);
    file_info(indx_del) = [];
end

% reads in image stack and represents it as a 3D array
for k = 1 : length(file_info)
    NUC(:,:,k) = imread(file_info(k).name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% read in CELL images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change current directory to CELL folder
cd ..
cd('./CELL')

% file extension
file_extension='*.png';
file_info=dir(file_extension);

% algorithm to delete pseudo-files
indx_del = [];
for k = 1 : length(file_info)
    filename=file_info(k).name;
    if startsWith(filename, '._')
        indx_del(k) = k;
    end
end

if exist('indx_del','var') == 1
    indx_del = nonzeros(indx_del);
    file_info(indx_del) = [];
end

% reads in image stack and represents it as a 3D array
for k = 1 : length(file_info)
    CELL(:,:,k) = imread(file_info(k).name);
end

%% 2) Segmentation of NUC
 
% change current directory back to parent directory 
cd ..

% Max intensity projection + contrast enhanced
MAX_NUC = imadjust(max(NUC,[],3));

% creates blue color scale for NUC
blues = linspace(0,1,max(unique(MAX_NUC)));
NUC_colors = zeros(max(unique(MAX_NUC)),3);
NUC_colors(:,3) = blues;

% applies 3D Gaussian filter to remove noise 
I2 = imgaussfilt3(NUC,2);

% Create nuclear segmentation

x = {}; % pre-allocation of variable for input statement
x{1} = 0; % pre-allocation of variable for input statement
nuc_thresh = 1; % preliminary threshold for NUC

while x{1} == 0

    % binarizes or segments the image stack
    nuc_bw = I2 > nuc_thresh;
    
    % plots original NUC image with binary mask
    figure
    subplot(1,2,1) % original
    imshow(MAX_NUC,NUC_colors); hold on;
    pbaspect([1 1 1])
    title('max intensity projection')
    subplot(1,2,2) % binary mask
    imshow(max(nuc_bw,[],3)); hold on;
    pbaspect([1 1 1])
    title(strcat('mask with threshold ={ }',num2str(nuc_thresh)))
    set(gcf, 'Position', get(0, 'Screensize'));

    % saves plot
    exportgraphics(gcf,strcat(filename(1:end-4),'_NUC_mask.png'),'Resolution',300)

    % handles input statements and iterative adjustments to threshold value
    x = inputdlg({'Segementation ok? Answer 1 for yes or 0 for no.','New Threshold Value. If you answered 1 previosly, type a random number.'},...
              'User Inputs', [1 75; 1 75]); hold off;

    close all

    % updates threshold value and x which controls if the while loop ends
    nuc_thresh = str2num(x{2});
    x{1} = str2num(x{1});

end

%% 3) Segmentation of CELL

% Max intensity projection + contrast enhanced
MAX_CELL = imadjust(max(CELL,[],3));

% creates red color scale for CELL
reds = linspace(0,1,max(unique(MAX_CELL)));
CELL_colors = zeros(max(unique(MAX_CELL)),3);
CELL_colors(:,1) = reds;

% applies 3D Gaussian filter to remove noise 
I2 = imgaussfilt3(CELL,2);

% Create cell segmentation
x = {}; % pre-allocation of variable for input statement
x{1} = 0; % pre-allocation of variable for input statement
CELL_thresh = 10; % preliminary threshold for CELL

while x{1} == 0

    % binarizes or segments the image stack
    CELL_bw = I2 > CELL_thresh;
    
    % plots original CELL image with binary mask
    figure
    subplot(1,2,1) % original
    imshow(MAX_CELL,CELL_colors); hold on;
    pbaspect([1 1 1])
    title('max intensity projection')
    subplot(1,2,2) % binary mask
    imshow(max(CELL_bw,[],3)); hold on;
    pbaspect([1 1 1])
    title(strcat('mask with threshold ={ }',num2str(CELL_thresh)))
    set(gcf, 'Position', get(0, 'Screensize'));

    % saves plot
    exportgraphics(gcf,strcat(filename(1:end-4),'_CELL_mask.png'),'Resolution',300)

    % handles input statements and iterative adjustments to threshold value
    x = inputdlg({'Segementation ok? Answer 1 for yes or 0 for no.','New Threshold Value. If you answered 1 previosly, type a random number.'},...
              'User Inputs', [1 75; 1 75]); hold off;

    close all

    % updates threshold value and x which controls if the while loop ends
    CELL_thresh = str2num(x{2});
    x{1} = str2num(x{1});

end

% creates a composite image
COMPOSITE(:,:,1) = MAX_CELL;
COMPOSITE(:,:,2) = zeros(size(MAX_CELL));
COMPOSITE(:,:,3) = MAX_NUC;


%% 4) Watershed segmentation

% Create watershed segmentation
x = {}; % pre-allocation of variable for input statement
x{1} = 0; % pre-allocation of variable for input statement
minima_cutoff = 4; % preliminary minima cutoff

while x{1} == 0

    % distance transform
    D = -bwdist(~CELL_bw);

    % converts all pixels that are not in CELL binary maks to Inf
    D(~CELL_bw) = Inf;

    % suppress regional minima in image using H-minima transform
    D = imhmin(D,minima_cutoff);

    % converts all pixels that are in NUC to -Inf
    D(nuc_bw) = -Inf;

    % perform watershed and set all background pixels to 0
    L = watershed(D);
    L (~CELL_bw) = 0;

    % pre-allocation for the creation of a random color scale
    colors = [0,0,0];
    for i = 1:max(unique(L))
        colors = [colors;rand(1,3)];
    end

    % plots original images with segmentation 
    figure (1)
    subplot(1,2,1) % original    
    imagesc(COMPOSITE); hold on;
    pbaspect([1 1 1])
    title('max intensity projection - contrast enhanced')
    subplot(1,2,2) % segmentation    
    imshow(max(L,[],3),colors); hold off;
    pbaspect([1 1 1])
    title(strcat('mask with minima cut-off ={ }',num2str(minima_cutoff)))
    set(gcf, 'Position', get(0, 'Screensize'));

    % saves plot
    exportgraphics(gcf,strcat(filename(1:end-4),'_CELL_segmentation_minima_thresh_',num2str(minima_cutoff),'.png'),'Resolution',300)

    % handles input statements and iterative adjustments to minima cut-off
    % value
    x = inputdlg({'Segementation ok? Answer 1 for yes or 0 for no.','New minima cut-off value. If you answered 1 previosly, type a random number.'},...
              'User Inputs', [1 75; 1 75]); hold off;

    close all

    % updates minima cut-off value and x which controls if the while loop ends
    minima_cutoff = str2num(x{2});
    x{1} = str2num(x{1});

end

% finds unique volumes in the NUC binary mask
CC = bwconncomp(nuc_bw);

% labels connected volumes
L_nuc = labelmatrix(CC);

% creates file names for saved NUC and CELL segmentation 
Lname = strcat(filename(1:end-4),'_CELL_SEGMENTATION.mat');
L_nuc_name = strcat(filename(1:end-4),'_NUC_SEGMENTATION.mat');

% save files
save(Lname,'L', '-v7.3')
save(L_nuc_name,'L_nuc', '-v7.3')

