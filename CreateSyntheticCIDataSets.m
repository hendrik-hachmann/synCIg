%% MIT License
% Copyright 2021 Hendrik Hachmann

% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so, 
% subject to the following conditions:

% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Create synthetic cochlear implant data-set
%
% Creates N volumes in an data structure
% There are several visualization and export functions 
% Per default CI data-sets and corresponding ground truth are exported
% as *.nii and *.nii.gz data

%% Parameters to control the properties of the data-set:
% A. Effect of region of interest
% B. Effect of voxel size
% C. Effect of left to right or right to left cochlea
% D. Effect of angles
% E. Effect of blobs
% F. Effect of number of electrodes 12, 16, 20, 22
% G. Effect of adding structures adding noise at pixel level
% H. Effect of intensity range

%% Create synthetic cochlear implant data set 
% reads ".nii.gz" files
% one folder should contrain modalities as 
% t1.nii.gz, flari.nii.gz, t1ce.nii.gz, t2.nii.gz
% and groundtruth data as: truth.nii.gz
       
% Create this number of volumes
dataSetSize = 5;

for dSS = 1:dataSetSize
    clearvars -except dSS
    fprintf('Creating dataset: %d\n', dSS);    
    
    %% Adjustable parameters
    % A. Effect of region of interest. 
    % ROI size as a factor of the CI dimensions
    guardIntervall = 10; % pixel

    % B. Effect of voxel siz
    % Voxel size ~ resolution -> gauss pyramid
    % Electrode position parameters
    distanceBetweenElectrodes = 4+randLH(-1,1)*2; % voxels -> determines volume resolution
    lengthVariation = 1; % factor
    angleVariation = 1; % factor

    % C. Effect of left to right or right to left cochlea
    mirror = true;
    theta_x = randLH(-1,1)*pi;
    theta_y = randLH(-1,1)*pi;
    theta_z = randLH(-1,1)*pi;

    % D. Effect of angles. 
    % Angles: 1) angle for helix in plane 2) deviation from plane (positiv z axis or negativ z axis)
    inPlaneAngleIncrease = 1.5+randLH(-1,1)*0.5;
    deviationFromPlane = 0.4+randLH(-1,1)*0.2;

    % E. Effect of blobs. 
    % Transition between nocontrast region vs. distantly space delectrode arrays
    % ~numberOf Smoothings: (imgausfilt3)
    sigma = 2+(rand-0.5); 
    noContrastRegionScore = 10;

    % F. Effect of number of electrodes 12, 16, 20, 22
    indexE = randi([1,4]);
    arrayE = [12, 16, 20, 22];
    numberOfElectrodes = arrayE(indexE);

    % G. Effect of adding structures adding noise at pixel level
    numberOfAdditionalBlobs = 10+randi([-9,9]);
    noisePower = 301+randi([-300,300]);
    minDistance2Electrodes = distanceBetweenElectrodes*1.5;

    % H. Effect of intensity range
    min2MaxIntesityValues = 4000+randi([-300,300]);

    %% Create data-set
    initPos1 = [1 ,1];         
    initPos2 = initPos1 + [distanceBetweenElectrodes,0];
    angel_pre = inPlaneAngleIncrease*pi/180 ;

    [contour, angel_inloop] = createContour(numberOfElectrodes, initPos1, initPos2, distanceBetweenElectrodes, angel_pre, lengthVariation, angleVariation);

    if mirror == true
        tmp = contour(:,1);
        contour(:,1) = contour(:,2);
        contour(:,2) = tmp;
    end

    z2 = linspace(1,numberOfElectrodes,numberOfElectrodes);
    z2 = z2 * deviationFromPlane;
    contour = [contour z2'];
    Matrix_rotation_x = [1 0 0;0 cos(theta_x) -sin(theta_x);0 sin(theta_x) cos(theta_x)];
    Matrix_rotation_y = [cos(theta_y) 0 sin(theta_y);0 1 0;-sin(theta_y) 0 cos(theta_y)];
    Matrix_rotation_z = [cos(theta_z) -sin(theta_z) 0; sin(theta_z) cos(theta_z) 0; 0 0 1];

    contour_rotation = contour * Matrix_rotation_z * Matrix_rotation_x * Matrix_rotation_y;
    contour_gt_tmp = contour_rotation;
    contourPel = round(contour_rotation);
    offset2GT = contourPel-contour_gt_tmp;

    % Spread points to blobs, blur successively 
    minX = min(contourPel(:,1))-1-guardIntervall;
    minY = min(contourPel(:,2))-1-guardIntervall;
    minZ = min(contourPel(:,3))-1-guardIntervall;
    contourPel2 = contourPel - [minX,minY,minZ];

    contourPel3 = contourPel2;
    contour_gt = contourPel3-offset2GT;

    maxX = max(contourPel3(:,1))+guardIntervall;
    maxY = max(contourPel3(:,2))+guardIntervall;
    maxZ = max(contourPel3(:,3))+guardIntervall;

    imageSize = [maxX maxY maxZ];
    I = zeros(imageSize);

    % Use upsampling and downsampling to get subvoxel positions
    upDownSamplingFactor = 5;
    I_up = imresize3(I, upDownSamplingFactor);
    contourPel_up = round(contour_gt*upDownSamplingFactor);
    linearInd = sub2ind(size(I_up), contourPel_up(:,1),contourPel_up(:,2),contourPel_up(:,3));
    I_up(linearInd) = 2^12;

    % Add false positives to data set
    if numberOfAdditionalBlobs > 0
        loop = true;
        falsePoints = [];
        while(loop)
            falsePointsCandidate = [randi(maxX*upDownSamplingFactor), randi(maxY*upDownSamplingFactor), randi(maxZ*upDownSamplingFactor)];

            dist2Electrodes = sqrt(sum((contourPel_up - falsePointsCandidate).^2,2));

            if min(dist2Electrodes) > minDistance2Electrodes*upDownSamplingFactor
                falsePoints = [falsePoints; falsePointsCandidate];
            end

            if size(falsePoints,1) >= numberOfAdditionalBlobs
                loop = false;
            end
        end
        linearInd2 = sub2ind(size(I_up), falsePoints(:,1),falsePoints(:,2),falsePoints(:,3));
        I_up(linearInd2) = 2^12;
    end

    for n = 1:noContrastRegionScore
            I_up(linearInd) = 2^12;
            if numberOfAdditionalBlobs > 0
                I_up(linearInd2) = 2^12;
            end
            I_up =  imgaussfilt3(I_up, sigma*upDownSamplingFactor,'padding',0);
    end

    % Downsample to original size
    I_down = imresize3(I_up,1/upDownSamplingFactor);

    I2 = I_down - min(I_down(:));
    currentRange = max(I2(:));
    spreadByFactor = min2MaxIntesityValues/currentRange;

    % Add noise
    [roiX, roiY, roiZ] = size(I2);
    additiveNoise = (rand(roiX, roiY, roiZ)-0.5)*2*noisePower;
    I3 = I2 ;
    contour_gauss = uint16(I3*spreadByFactor + additiveNoise);


    %% Visualize data
    if 0
        % Show CI in 2D
        figure();
        plot(contour(:,1), contour(:,2),'r');
        hold on;       
        plot(contour(:,1), contour(:,2),'r*');  
    end

    if 0
        % Show CI in 3D
        figure();
        plot3(contour_rotation(:,1),contour_rotation(:,2),contour_rotation(:,3), 'r*');
        hold on;
        plot3(contour_rotation(:,1),contour_rotation(:,2),contour_rotation(:,3), 'r');
        axis equal
    end

    if 0
        % Show CI cropped and  with guardintervall in 3D
        figure();
        plot3(contourPel3(:,1),contourPel3(:,2),contourPel3(:,3), 'r*');
        hold on;
        plot3(contourPel3(:,1),contourPel3(:,2),contourPel3(:,3), 'r');
        if numberOfAdditionalBlobs > 0
            plot3(falsePoints(:,1), falsePoints(:,2), falsePoints(:,3), 'b*');
        end
    end

    if 0
        % Show CI slice by slice
        cLim = [1, 2^12];
        for nn = 1:size(contour_gauss,3)    
            figure(12);
            hold off;
            imagesc(contour_gauss(:,:,nn),cLim);     
            hold on;
            imagesc(contour_gauss(:,:,nn),cLim);     
            colorbar;
            pause(0.05);
        end
    end

    %% Export data
    if 1
        if 0
            % Save slice through as video    
            folderNameTmp = './tmp/';    
            listing = dir([folderNameTmp]);
            n1 = size(listing,1) + 1;
            [status, msg, msgID] = mkdir([folderNameTmp,sprintf('/vid%d/',n1)]);
            cLim = [1, 2^12];
            for nn = 1:size(contour_gauss,3)  
                h = figure('visible','off');
                imagesc(contour_gauss(:,:,nn),cLim);
                colorbar;
                set(gcf, 'paperpositionmode', 'auto'); % to prevent matlab crashing
                filename_h = [folderNameTmp, sprintf('vid%d/%04d',n1,nn), '.png'];
                saveas(h,filename_h);      
            end
        end

        % Make video
        if 0
            command = sprintf('cat $(ls %s/*.png) | ffmpeg2 -r 25 -i - -c:v libx264 -vf scale=1920:1080 -pix_fmt yuv420p %s/sliceVolume.mp4 -y', [folderNameTmp,sprintf('/vid%d/',n1)], [folderNameTmp,sprintf('/vid%d/',n1)]);
            system(command);        
            command2 = sprintf("echo '#!/bin/bash'  > %sstartSliceVolume.sh", [folderNameTmp,sprintf('/vid%d/',n1)]);
            system(command2);      
            command3 = sprintf("echo 'mplayer -fs -fps 5 sliceVolume.mp4'  >> %sstartSliceVolume.sh", [folderNameTmp,sprintf('/vid%d/',n1)]);
            system(command3);      
            command4 = sprintf('chmod 770 %sstartSliceVolume.sh', [folderNameTmp,sprintf('/vid%d/',n1)]);
            system(command4);      
            command5 = sprintf("echo 'xdotool key ctrl+shift+q' >> %sstartSliceVolume.sh", [folderNameTmp,sprintf('/vid%d/',n1)]);
            system(command5);      
        end

        % Save volume and ground truth as *.mat
        if 0
            dataSet = contour_gauss;
            save([folderNameTmp,sprintf('/vid%d/',n1), 'volume.mat'], 'dataSet');
            groundTruth = contour_gt;
            save([folderNameTmp,sprintf('/vid%d/',n1),'groundtruth.mat'], 'groundTruth');
        end

        % Export data set as .nii.gz
        if 1           
            folderNameTmp = './tmp/dataset/';    
            currendDataFolder = [folderNameTmp, sprintf('/dataset%d/',dSS)];
            [status, msg, msgID] = mkdir([folderNameTmp, sprintf('/dataset%d/',dSS)]);
       
            dataSet =  imresize3(contour_gauss, [114 114 114]);
            
            niftiwrite(dataSet, [currendDataFolder,'/t1.nii']);     
            command = ['gzip -c ', currendDataFolder, '/t1.nii > ', currendDataFolder, '/t1.nii.gz'];
            system(command);
            
            groundTruth = contour_gt;
            gtVolume = 0 * contour_gauss;
            origSize = size(gtVolume);
            targetSize = [114 114 114];
            gtVolume =  imresize3(gtVolume, targetSize);   
            
            transformationFactor =  targetSize./origSize;
            groundTruth_TargetSize = groundTruth .* transformationFactor;
            idx = sub2ind(size(gtVolume), round(groundTruth_TargetSize(:,1)),round(groundTruth_TargetSize(:,2)),round(groundTruth_TargetSize(:,3)));                          
            gtVolume(idx) = 1;       
            
            se = strel('sphere',2);
            gtVolume2 =  imdilate(gtVolume,se);
            
            niftiwrite(gtVolume2, [currendDataFolder,'/truth.nii']);     
            command = ['gzip -c ', currendDataFolder, '/truth.nii > ', currendDataFolder, '/truth.nii.gz'];
            system(command);
        end
    end
end

disp('End');

%% Create Contour points
function [contour,angel_inloop] = createContour(n,  initPos1, initPos2, distance2Neighboor, angel_pre, len_variation, angel_variation)
    contour(1,:) = initPos1;
    contour(2,:) = initPos2;
    angel_inloop = 0;
    for i =3:n
        angel_inloop = angel_inloop + angel_pre * i * angel_variation;       
        contour(i,1) = contour(i-1,1) + cos(angel_inloop) * distance2Neighboor * len_variation;
        contour(i,2) = contour(i-1,2) + sin(angel_inloop) * distance2Neighboor * len_variation;
    end
end

%% Random sample from interval
function [outputArg1] = randLH(inputArg1, inputArg2)
    diffLH =  inputArg2 - inputArg1;
    tmp = rand*diffLH;
    outputArg1 = tmp-diffLH/2;
end
