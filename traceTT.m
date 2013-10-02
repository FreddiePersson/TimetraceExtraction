function [tracedTT TTspecs] = traceTT(TT, varargin)
% function [tracedTT TTspecs] = alignTT(TT)
%   Assumes center of molecule approximatly in the center of timetrace and
%   with a 10 pixels space on each edge. TT is a matrix that can contain multiple 
%   timetraces (of the same size) so that, TT(:, :, i) is the i:th timetrace.

% If nargin > 1 the average of all timetraces should be traced but all
% aligned individually. Else each is treated individually.
if nargin>1
    singleTrace = true;
    numTT = 1;
else
    singleTrace = false;
    numTT = size(TT, 3);
end

% Tolerance of the length in relation to the previos frame/line
pulseTol = 0.8;
% Number of pulses
nPulses = 1;

% Length of fitting in the channel. Normally the whole span.

% The initial guess of the pulse height related to the max value of the
% first frame/line in each TT
pulseHeight = 0.9;


% Common parameters for the TTs
chanL = size(TT, 2);
x = 1:chanL;
nFrames = size(TT, 1);

% Initial approximation of edges etc
temp = mean(TT, 3);
temp = temp(1, :);
f = figure;
plot(temp, '-r')
[px ~] = ginput(2)
close(f); clear temp;
px1 = min(px);
px2 = max(px);
x0 = (px1 + px2)/2;
L0 =  px2 - px1; 
sigma = 20;
clear px px1 px2;

% Length of fitting in the channel. Normally the whole span.
sx1 = 1;
sx2 = chanL;

% Initiate variables
currentTT = zeros(size(TT(:, :, 1)));
tracedTT = zeros(size(TT));
TTspecs = zeros(size(TT, 1), 2, numTT);  % xCent, Lx, TT

pause on;
for nTT = 1:numTT
    
    if singleTrace
        currentTT = mean(TT, 3);
    else
        currentTT = TT(:, :, nTT)
    end
    
    for frame = 1: nFrames
        intensity = currentTT(frame, :);
        intensity = (intensity-min(intensity))./max(intensity);
        if frame == 1

            % Set initial pulse params.
            p(1:2,1) = 0;	% I1 & I2
            p(1,2) = x0;	% x0
            p(2,2) = sigma;	% sigma
            p(1,3:3+nPulses-1) = pulseHeight;% A0
            p(2,3:3+nPulses-1) = L0;%L0;
            
            % Distribute pulses if nPulses > 1
            if (nPulses > 1)
                % Lx
                p(2,3:3+nPulses-1) = (p(2,3:3+nPulses-1) ./ sum(p(2,3:3+nPulses-1))) .* L0;
                % Ax
                p(1,4:2:2+nPulses) = .6*pulseHeight;% A0
            end
            
            % Initialize result data
            prevFrame = 0;
            chk = true;
        end
        
        % Fit parameters to data
		warning off;
		p = nlinfit(x,intensity,@modelPulse,p);
		warning on;
		
		if (frame == 1)
			pFirst = p;
        end
        
        % Get fitted pulse
		fittedPulse = modelPulse(p,x);
        plot(intensity, '-r')
        hold on
        plot(fittedPulse);
        hold off
%         print(gcf, '-dpdf', num2str(frame));
        pause(0.05)
        
		% Get molecule meta data
		TTspecs(frame, 1, nTT) = p(1,2);
		TTspecs(frame, 2, nTT) = sum(p(2,3:end));
		
        % Check molecule within tolerance
		if (prevFrame > 0)
			if ( abs(TTspecs(frame, 2, nTT)-TTspecs(prevFrame, 2, nTT)) > pulseTol*TTspecs(prevFrame, 2, nTT))
				warning(['Trace ' num2str(nTT) ' of ' num2str(size(TT, 3)) ' failed due to length variation between frames are out of tolerance'])
				break;
			end
        end
        
		% Counter
		prevFrame = frame;
		
	end
        

end
pause off;

for nTT = 1:size(TT, 3)
    if singleTrace
        temp = alignTT(TT(:, :, nTT), TTspecs(:, :, 1)); 
    else
        temp = alignTT(TT(:, :, nTT), TTspecs(:, :, nTT));
    end
    tracedTT(:, :, nTT) = temp;   
end
close(gcf);
end






function alignedTT = alignTT(TT, specs)
% Aligns a single timetrace (TT).

% Molecule center
xCent = specs(:, 1);

% Chooses suitable padding value
pad = min(min(TT));

shiftPos = round(xCent-xCent(1));
alignedTT = TT;
for i = 1:size(TT, 1)
    
	temp = TT(i, :);
    if shiftPos(i) == 0
        
    elseif shiftPos(i)>0
        temp(1:abs(shiftPos(i))) = [];
        temp = [temp, ones(1, abs(shiftPos(i)))*pad]; 
        alignedTT(i, :) = temp;
    else
        temp(end-abs(shiftPos(i))+1:end) = [];
        temp = [ones(1, abs(shiftPos(i)))*pad, temp];
        alignedTT(i, :) = temp;
    end
end


% 	% Write to tiff file
% 	imwrite(...
% 		uint16(alignedTT),...
% 		[ savePathName, strtok(saveFileName,'.') '.tif'],...
% 		'WriteMode','overwrite');
	
end


