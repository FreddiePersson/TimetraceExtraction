function alignedStack = alignStack(stack, specs)
% Aligns a single timetrace (TT).

% Molecule center
xCent = specs(:, 1);

% Chooses suitable padding value
pad = min(min(stack(:, :, 1)));

shiftPos = round(xCent-xCent(1));
alignedStack = stack;
for i = 1:size(stack, 3)
    
    for ii = 1:size(stack, 1)
        temp = stack(ii, :, i);
        if shiftPos(i) == 0
            
        elseif shiftPos(i)>0
            temp(1:abs(shiftPos(i))) = [];
            temp = [temp, ones(1, abs(shiftPos(i)))*pad];
            alignedStack(ii, :, i) = temp;
        else
            temp(end-abs(shiftPos(i))+1:end) = [];
            temp = [ones(1, abs(shiftPos(i)))*pad, temp];
            alignedStack(ii, :, i) = temp;
        end
    end
%     i 
%     shiftPos(i)
end


% 	% Write to tiff file
% 	imwrite(...
% 		uint16(alignedTT),...
% 		[ savePathName, strtok(saveFileName,'.') '.tif'],...
% 		'WriteMode','overwrite');
	
end
