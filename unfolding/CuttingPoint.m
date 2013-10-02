function [] = CuttingPoint(fileName)

    
    % Load the image from .tif file.
    imgArr = GetImgArrFromFile(fileName);    
    imgArr = imgArr ./ max(imgArr(:));
    origImgArr = imgArr;
    
    
    % Smooth the image to improve region contrast.
    imgArr = smoothimg(imgArr);
    
    
    % Thresholding via Otsu's method. seg_I will be a binary image mask for
    % the segmented image.
    thresh = multithresh(imgArr,2);
    seg_I = imquantize(imgArr, thresh);
    linearDNAMask = (seg_I == 2);
    linearDNAMask = imfill(linearDNAMask,'holes');
    circleDNAMask = (seg_I == 3);
    circleDNAMask = imfill(circleDNAMask,'holes');
    seg_I(seg_I~=3) = 0;
    seg_I(seg_I==3) = 1;
    
    
    % Do the thresholding again if we incorrectly identified the entire DNA
    % region as the "light" region. This sometimes happens if there's too 
    % much background noise. We say that this scenario happened if the
    % bottom k rows all contain "light region."
    k = 10; 
    if sum(zeros(k+1,1) < sum(seg_I(end-k:end,:),2)) == (k + 1)
        imgArr(~seg_I) = 0;
        thresh = multithresh(imgArr,2);
        seg_I = imquantize(imgArr, thresh);
        linearDNAMask = (seg_I == 2);
        linearDNAMask = imfill(linearDNAMask,'holes');
        circleDNAMask = (seg_I == 3);
        circleDNAMask = imfill(circleDNAMask,'holes');
        seg_I(seg_I~=3) = 0;
        seg_I(seg_I==3) = 1;
    end
    
    
    % Destroy everything except the largest connected component.
    CC = bwconncomp(seg_I);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    seg_I(CC.PixelIdxList{idx}) = 2;
    seg_I(seg_I ~= 2) = 0;
    seg_I(seg_I == 2) = 1;
    
    
    % Trim edges of the main region. leftPoint and rghtPoint are the left
    % and right columns marking the boundaries of the circular DNA region.
    % regionLoc is the rough guess for the breaking time.
    [seg_I, leftPoint, rghtPoint, regionLoc] = TrimEdges(seg_I);
        
     
    % Find the breaking point by finding the tallest black bar under the
    % top light region.
    heights = zeros(rghtPoint-leftPoint,1);
    for i=1:length(heights)
        lastBrightPt = find(seg_I(:,i+leftPoint),1,'last');
        if ~isempty(lastBrightPt)
            heights(i) = sum(seg_I(regionLoc:end,i+leftPoint) == 0);
        end
        heights(heights == 0) = [];        
    end
    heightsBreakIdx = round(mean(find(heights == max(heights))))-1;
    breakPrctVal = heightsBreakIdx / (length(heights)-1);
    
    breakXCoord = heightsBreakIdx+leftPoint+1;
    breakYCoord = find(seg_I(regionLoc:end,breakXCoord) == 1, 1, 'last')+regionLoc;
        
    
    % Now get the mean widths of the linear and circular DNA regions
    circularRegionBottom = find(sum(seg_I,2) == 0, 1, 'first');
    meanLinearLength = mean(sum(linearDNAMask(circularRegionBottom:end,:),2));
    meanCircleLength = mean(sum(circleDNAMask(1:breakYCoord,:),2));
    
    
    % Finally get the "unfolding time," defined as the number of rows
    % between the breaking point and the end of the circular region.
    unfoldingTime = circularRegionBottom - breakYCoord;
    
    
    % Save the image.
    DispMontage(origImgArr, seg_I, breakXCoord, breakYCoord);
    
    
    % Save the output data in a struct.
    OutputData = struct();
    OutputData.BreakPosition = breakPrctVal;
    OutputData.BreakTime = breakYCoord;
    OutputData.CircleDNALength = meanCircleLength;
    OutputData.LinearDNALength = meanLinearLength;
    OutputData.UnfoldingTime = unfoldingTime;
    
    % Pring a message.
    disp(['Break position:' 9 9 num2str(OutputData.BreakPosition)])
    disp(['Break time:' 9  9 num2str(OutputData.BreakTime)])
    disp(['Circle DNA Length:' 9 num2str(OutputData.CircleDNALength)])
    disp(['Linear DNA Length:' 9 num2str(OutputData.LinearDNALength)])
    disp(['Unfolding time:' 9 9 num2str(OutputData.UnfoldingTime)])
    
    savefname = strcat(fileName(1:end-4), '_OUTPUT.mat');
    save(savefname,'OutputData');    
    
    disp(['Data saved to:' 9 9 savefname])
    

end


%% Subfunctions
function [ imgArr ] = GetImgArrFromFile(fileName)

    tiffInfo = imfinfo(fileName);   % Get the TIF file information
    no_frame = numel(tiffInfo);     % Get the number of images in the file
    traceCell = cell(no_frame,1);   % Preallocate the cell array
    for frame = 1:no_frame
      traceCell{frame} = im2double(imread(fileName,'Index',frame,'Info',tiffInfo));
    end
    imgArr = traceCell{1,1};
    
end


function [ outArr ] = smoothimg(imgArr)

% Use a 2D moving average because this is the optimal smoothing technique
% for uncovering trends obscured by normally-distributed fluctuations.
f = fspecial('average',[7 3]);
outArr = conv2(imgArr, f, 'valid');

end


function [ seg_I, lp, rp, regionLoc ] = TrimEdges(seg_I)
   
    
    % Find the vertical coordinate where the breaking region begins.
    % Consider the horizontal sum. I.e., the width of the light region for
    % a given row. Then it'll begin high and stay flat until the breaking
    % region where it drops very quickly (or more slowly, depending on the
    % speed of the break). Then it remains zero until the bottom of the
    % image. What we'll do is define this sum and eliminate the
    % zero-region. Then at each point, we'll compute the variance of all
    % following points. This value should reach its maximum just before the
    % drop-off.
    s = sum(seg_I,2);
    s(s == 0) = [];
    v = zeros(size(s));
    for i = 1:length(v)
        v(i) = var(s(i:end));
    end
    lcInit = find(v == max(v),1,'first');    
    lc = lcInit - floor(size(seg_I,1) / 10);
    if lc < 1
        lc = 1;
    end
    
    
    % Find the left and right points of the light region. These will be the
    % farthest right and farthest left (respectively) black points found
    % within the light region above lc, the vertical coordinate that we
    % identified as the beginning of the breaking region.
    leftPoint = 1;
    rghtPoint = inf;
    for i = 1:lc
        lp = find(seg_I(i,:),1,'first');
        rp = find(seg_I(i,:),1,'last');
        if lp > leftPoint
            leftPoint = lp;
        end
        if rp < rghtPoint
            rghtPoint = rp;
        end
    end
    
    
    % Destroy all the light points outside [lp, rp].
    seg_I(:,1:lp) = 0;
    seg_I(:,rp:end) = 0;
    
    
    % Fill holes.
    seg_I = imfill(seg_I,'holes');
    se = strel('line',3,90);
    seg_I = imclose(seg_I,se);
    
    % Will return lp and rp as-is, so we just need to define output
    % argument regionLoc before ending the function.
    regionLoc = lc;
    
end


function [] = DispMontage(origImgArr, seg_I, breakXCoord, breakYCoord)

% Generate the figure.
figure();

% Put a circle on top of the breaking point in the image.
cx=breakXCoord;
cy=breakYCoord;
ix=size(origImgArr,2);
iy=size(origImgArr,1);
r=2;
[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((x.^2+y.^2)<=r^2);

% Make the circle red.
origImgArr = repmat(origImgArr, [1 1 3]);
origImgArrRed = origImgArr(:,:,1);
origImgArrGrn = origImgArr(:,:,2);
origImgArrBlu = origImgArr(:,:,3);
origImgArrRed(c_mask) = 1;
origImgArrGrn(c_mask) = 0;
origImgArrBlu(c_mask) = 0;
origImgArr(:,:,1) = origImgArrRed;
origImgArr(:,:,2) = origImgArrGrn;
origImgArr(:,:,3) = origImgArrBlu;

% Generate the side-by-side montage with the original kymograph on the left
% and the mask on the right.
imshowpair(origImgArr, seg_I, 'montage');

end