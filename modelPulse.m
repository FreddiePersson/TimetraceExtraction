function y = modelPulse(p,x,varargin)
	% MODELPULSE Creates a pulse with variable number of features.
	%
	% Syntax:
	%
	%   y = modelPulse(p,x)
	%   y = modelPulse(p,x,'debug')
	%
	% Description:
	%
	%   modelPulse(p,x) Creates a pulse on x with variable 
	%   number of features n defined by p.
	%   p is defined as:
	%   
	%     p = 
	%         
	%         I1  x0  [Ax(1) Ax(2) ... Ax(n)]
	%         I2  sx  [Lx(1) Lx(2) ... Lx(n)]
	%   
	%   Where:
	%     
	%     I1, I2   Initial and ending intensity value before/after
	%              pulse.
	%     x0       Center of pulse on x.
	%     sx       Slope of all features in pulse.
	%     Ax       Amplitude of features
	%     Lx       Width of features. Total width of pulse is 
	%              constrained to within x.
	%   
	%   modelPulse(...,'debug') prints variables as a struct in 
	%   prompt.
	%   
	
	function xHat = reOdd(x,L)
		
		s = (-1).^(floor(x./L));
		xHat = s .* (mod(x,L) - L.*(1 - s)./2 );
		
	end
	
	% Get and calc. parameters
	dx = length(x);
	nFeatures = size(p,2) - 2;
	
	I1 = reOdd(p(1,1),1); % Within unity
	I2 = reOdd(p(2,1),1); % Within unity
	
	x0 = reOdd(p(1,2),dx);				% pulse center within dx
	minD = min(max(x)-x0,x0-min(x));	% Shortest dist. to edge
	
	sx = dx ./ (10*(3+reOdd(p(2,2),10)));	% Scale with dx
	
	Ax = [I1 reOdd(p(1,3:end),1) I2]; 		% Within unity
	
	Lx = [0 p(2,3:end)];		% Add zero pulse width first
	if (sum(Lx)/2 > minD)	% Pulse width within dx
		Lx = (Lx ./ sum(Lx)) .* minD;
	end
	
	
	% Calc. function values
	L0 = sum(Lx);
	z = ( x - (x0 - L0/2) ) ./ ( sqrt(2) .* sx );
	y = zeros(size(x))+I1;
	Lx = cumsum(Lx);
	z0 = Lx ./ ( sqrt(2) .* sx );
	a = (Ax(2:end) - Ax(1:end-1)) ./ 2;
	
	% Add pulses
	for i = 1 : nFeatures+1
	
		y = y + a(i).*(1 + erf(z - z0(i)));
		 
	end

	% Debug display
	if (nargin == 3)
		if (strcmp(varargin{1},'debug'))
		
			m.dx = dx;
			m.nFeatures = nFeatures;
			m.I1 = I1;
			m.I2 = I2;
			m.x0 = x0;
			m.minD = minD;
			m.L0 = L0;
			m.Lx = Lx;
			m.Ax = Ax;
			m.sx = sx;
			m.a = a;
			
			disp(m)
		end
	end
	
end