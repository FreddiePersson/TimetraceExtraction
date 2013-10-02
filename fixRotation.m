function [rotData, rotAngle] = fixRotation(rawData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


maxFrame = mean(rawData, 3);

figure('Name', 'Choose a region');imagesc(maxFrame);
rect = getrect(gcf);
close(gcf)
x1 = ceil(rect(1)); x2 = floor(x1+rect(3)); y1 = ceil(rect(2)); y2 = floor(y1+rect(4));
maxFrameCrop = maxFrame(y1:y2, x1:x2);

maxFrameCrop = mat2gray(maxFrameCrop);
levelOtsu = graythresh(maxFrameCrop);
BW = im2bw(maxFrameCrop, levelOtsu);
molPeaks = BW.*maxFrameCrop;
% figure; imagesc(molPeaks);

[yData xData] = find(BW);

warning off
rrFit = robustfit(xData, yData);
warning on

rotAngle = atand(rrFit(2))

rotData = imrotate(rawData, rotAngle, 'bicubic', 'crop');

rotData(find(rotData < min(rawData(:)))) = min(rawData(:));

figure;imagesc(mean(rotData, 3));
disp('Finito')
 
	
end

