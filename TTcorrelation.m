function [lineCorrelations, alignedData] = TTcorrelation(varargin)

%% Function to check correlation between top and bottom halves of TTs.
% For single row TTs place the different rows in a tiff stack and run.
% This is a test so that the result can be modified after what is found to
% be best. The general strategy is at least to compare correlations...

display = false
chooseCorrROI = false
ChannelWidth = 6  % Put to 0 if dont want to use it...

%% Options
timeAverage = 5;
if length(varargin)>0 & isnumeric(varargin{1})
timeAverage = varargin{1};
end

%% Code

loadStack;

% Fix rotation for the rawData stack
[rotData, rotAngle] = fixRotation(rawData);

% Gives aligned but not stretched TTs back... If stack of TTs is sent in
% then the sum of those is used for the alignment.
[tracedTTs, TTspecs] = extractLineTimetraces(rotData, 0, 6);


alignedData = alignStack(rawData, TTspecs);

temp = sum(alignedData, 3)./size(alignedData, 3);

if display
    figure; imagesc(temp);
end

if timeAverage >0
    lineCorrelations = zeros(size(tracedTTs, 1)/timeAverage, size(tracedTTs, 3), size(tracedTTs, 1)/timeAverage);
else
    lineCorrelations = zeros(size(tracedTTs, 1), size(tracedTTs, 3), size(tracedTTs, 1));
end

leftPix = round(TTspecs(1, 1))-ceil(TTspecs(1, 2)/2);
rightPix = round(TTspecs(1, 1))+ceil(TTspecs(1, 2)/2);
for ind = 1:size(tracedTTs, 3)
    temp = tracedTTs(:, :, ind);
    if timeAverage >0
        % Fix means along the time direction in the TT
        % p*q pieces of m*n blocks
        m = timeAverage;
        n = 1;
        p = size(temp, 1)/m;
        q = size(temp, 2)/n;
        
        [I,J,K,L] = ndgrid(1:m,0:n-1,0:p-1,0:q-1);
        temp = reshape(mean(reshape(temp(I+m*p*J+m*K+m*p*n*L),m*n,p*q),1),p,q);
    end
    
    for circInd = 1:size(temp, 1)
        % shift the matrix so row 2 becomes row 1 etc...
        temp = circshift(temp, -(circInd-1));
    for ind2 = 1:size(temp, 1)
        template = temp(1, leftPix:rightPix)-mean(temp(1, leftPix:rightPix));
        lineData = temp(ind2, leftPix:rightPix)-mean(temp(ind2, leftPix:rightPix));
       lineCorrelations(ind2, ind, circInd) = sum(template.*lineData)/(std(template).*std(lineData).*length(template));
    end
   end
end

if display
    % display the average correlations
    figure('Name', 'Mean of possible cross correlations');
    avLineCorr = nanmean(lineCorrelations, 3);
    surf(avLineCorr);
end

% remove 'autocorrelation', correlation between the same row and itself.
lineCorrelations(1, :, :) = [];
avLineCorr(1, :, :) = [];

if chooseCorrROI
figure('Name', 'Choose region to calculate mean and std'); imagesc(avLineCorr);
rect = getrect(gcf);
close(gcf);
% Check so that we dont end up with 0 or to large values
if round(rect(1))==0
    rect(1) = 1;
end
if round(rect(2))==0
    rect(2) = 1;
end
if round(rect(1)+rect(3))>size(avLineCorr, 1)
    rect(3) = rect(3)-1;
end
if round(rect(2)+rect(4))>size(avLineCorr, 2)
    rect(4) = rect(4)-1;
end

lineCorrCrop = lineCorrelations(round(rect(2)):round(rect(2)+rect(4)), round(rect(1)):round(rect(1)+rect(3)), :);

end

% Throw important things to the base workspace starting with the filename
[~, fnOrig, ~] = fileparts(full_filename);
fn = genvarname(fnOrig);
assignin('base', ['LCmat_' fn], lineCorrelations);
assignin('base', ['lineTTspecs_' fn], TTspecs);
assignin('base', ['lineTTs_' fn], tracedTTs);
assignin('base', ['rot_' fn], rotAngle);

disp(['Mean of correlation of total area: ' num2str(nanmean(lineCorrelations(:)))]);
disp(['Std of correlation of total area: ' num2str(nanstd(lineCorrelations(:)))]);

if chooseCorrROI
    disp(['Mean of correlation of cropped area: ' num2str(nanmean(lineCorrCrop(:)))]);
    disp(['Std of correlation of cropped area: ' num2str(nanstd(lineCorrCrop(:)))]);
    assignin('base', ['croppedLCmat_' fn], lineCorrCrop);
    assignin('base', ['mean_' fn], nanmean(lineCorrCrop(:)));
    assignin('base', ['std_' fn], nanstd(lineCorrCrop(:)));
    
else
    assignin('base', ['mean_' fn], nanmean(lineCorrelations(:)));
    assignin('base', ['std_' fn], nanstd(lineCorrelations(:)));
end
    

% save('temp.mat')
disp('Finito')

end
%% 

function alignedStack = alignStack(rawData, specs)
% Aligns a stack of images.

% Molecule center
xCent = specs(:, 1);

% Chooses suitable padding value
pad = min(min(min(rawData)));

shiftPos = round(xCent-xCent(1));
alignedStack = rawData;
for i = 1:size(rawData, 3)
    
	temp = rawData(:, :, i);
    if shiftPos(i) == 0
        
    elseif shiftPos(i)>0
        temp(:, 1:abs(shiftPos(i))) = [];
        temp = [temp, ones(size(rawData, 1), abs(shiftPos(i)))*pad]; 
        alignedStack(:, :, i) = temp;
    else
        temp(:, end-abs(shiftPos(i))+1:end) = [];
        temp = [ones(size(rawData, 1), abs(shiftPos(i)))*pad, temp];
        alignedStack(:, :, i) = temp;
    end
end
end

