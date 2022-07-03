function [amostra, sampleIndex] = sampleGaussianInterval(I, g, numSamples)
	
	amostra = nan(1,numSamples);
	if(g.s < 1e-6)
		sampleIndex = 1;
	else
		sampleIndex = 0;
		numIntervals = size(I,1);
		IntervalG = zeros(numIntervals,2);
		for k = 1 : numIntervals
			x = I(k,:);
			IntervalG(k,:) = g.A*exp(-0.5*((x - g.mu)/g.s).^2);
		end
		gII = [IntervalG I];
		[gII,~,~] = unique(gII, 'rows');
		If = gII(end,3:4);
		
		A = linspace(If(1), If(2), numSamples);
		gA = g.A*exp(-0.5*((A - g.mu)/g.s).^2);
		A = A';
		gA = gA';
		gAA = [gA A];
		[gAA,~,~] = unique(gAA, 'rows');
		amostra = gAA(numSamples : -1 : 1, 2);
	end
end
