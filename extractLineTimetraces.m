% Function to extract line timetraces in subsequent pixelrows 
function [alignedTTs, TTspecs] = extractLineTimetraces(rawData, rotAngle, channelWidth)
% rotAngle = 1;
align = true;
display = false;
tiltBaseline = false;
%channelWidth = 0; %In pixels. If set to 0 a ROI will be requested, else the
                    %channelWidths bightest pixelrows be used

if isempty(rawData)
loadStack;
end

if rotAngle~= 0
rotData = imrotate(rawData, rotAngle, 'bicubic', 'crop');
else
    rotData=rawData;
end

fhand = figure;
imagesc(rotData(:, :, 1));

if tiltBaseline
    % Draw a larger area than you want to extract TT from
    rect = getrect(fhand);
    rotDataCrop = rotData(floor(rect(2)):ceil(rect(2)+rect(4)), floor(rect(1)):ceil(rect(1)+rect(3)), :);%zeros(rect(3), rect(4), size(rotData, 3));
    
   %J = roifill();
   % Crop and flatten the data
   for i=1:size(rotDataCrop, 3)
       J=roifill(rotDataCrop(:, :, i), [1, 1, size(rotDataCrop, 2), size(rotDataCrop, 2)], [1, size(rotDataCrop, 1), size(rotDataCrop, 1), 1]);
       rotDataCrop(:, :, i)=rotDataCrop(:, :, i)-J;
   end
   rotData = rotDataCrop;
   fhand = figure;
   imagesc(rotData(:, :, 1));
end

if channelWidth == 0
    rect = getrect(fhand);
    
    % Check so that we dont end up with 0 or to large values
    if round(rect(1))==0
        rect(1) = 1;
    end
    if round(rect(2))==0
        rect(2) = 1;
    end
    if round(rect(1)+rect(3))>size(rotData, 1)
        rect(3) = rect(3)-1;
    end
    if round(rect(2)+rect(4))>size(rotData, 2)
        rect(4) = rect(4)-1;
    end
    
    x1 = round(rect(1)); x2 = round(x1+rect(3)); y1 = round(rect(2)); y2 = round(y1+rect(4));
    
else
    
    x1 = 1; x2 = size(rotData, 2);
    
    % Find the 'channelWidth' brightest pixel rows
    zProj = mean(rotData, 3);
    sumRows = sum(zProj, 2);
    maxVal = 0; maxInd = 0;

    for ind = 1:length(sumRows)-(channelWidth-1)
        if maxVal < sum(sumRows(ind:ind+channelWidth-1))
            maxVal = sum(sumRows(ind:ind+channelWidth-1));
            maxInd = ind;
        end
    end
    y1 = maxInd;
    y2 = y1+channelWidth-1;
    if display
        figure;hold on;plot(sumRows);
        plot([y1, y2], [max(sumRows) max(sumRows)], '+r'); hold off;
    end
end
close(fhand)

lineTT = zeros(1, (x2-x1)+1, size(rotData, 3));
rawTTs = zeros(size(rotData, 3), (x2-x1)+1, (y2-y1));
    
for ind = y1:y2
    lineTT = rotData(ind, x1:x2, :);
    lineTT = reshape(lineTT, [ceil(x2-x1)+1, size(rotData, 3)]);

    rawTTs (:, :, ind-y1+1) = lineTT';

end

if align
    [alignedTTs TTspecs] = traceTT(rawTTs, 'avTT');
end

if display
    
    if align
        figure('Name', ['Aligned mean TT']); imagesc(mean(alignedTTs, 3));
        for ind = 1:size(alignedTTs, 3)
            figure('Name', ['Aligned TT row ' num2str(ind) ' of ' num2str(size(alignedTTs, 3))]);
            imagesc(alignedTTs(:, :, ind));
        end
    else
        figure('Name', ['Raw mean TT']); imagesc(mean(rawTTs, 3));
        for ind = 1:size(rawTTs, 3)
            figure('Name', ['Raw TT row ' num2str(ind) ' of ' num2str(size(alignedTTs, 3))]);
            imagesc(rawTTs(:, :, ind));
        end
    end
end

% close all

% 
%  v = strrep(full_filename, '.tif', '_extension.csv');
%     dlmwrite(v,TTspecs);
%     w = strrep(full_filename, '.tif', '_traces.csv');
%     dlmwrite(w,alignedTTs);
% 
% 
% 	% Write to tiff file
% 	imwrite(...
% 		uint16(mean(rawTTs, 3)),...
% 		[ '/Users/FredrikPersson/Downloads/TTfromXXX.tif'],...
% 		'WriteMode','overwrite');
end