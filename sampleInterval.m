function [branches, branchNum] = sampleInterval(I, G, isArc, numSamples, myZero)
	
	I = juntaIntervalos(I, myZero);
	Ii = Interval2Struct(I, myZero);
	% oderna as gaussianas e escolhe as 2 de maior amplitude
	numGauss = length(G);
	vecMu = zeros(numGauss,2);
	for j = 1 : numGauss
		vecMu(j,:) = [G(j).A, j];
	end
	[vecMu,~,~] = unique(vecMu, 'rows');
	clear Gaux
	for j = 1 : min(numGauss,2)
		Gaux(j) = G(vecMu(numGauss-(j-1),2));
	end
	G = Gaux;

	numIntervals = length(Ii);
	branches = nan(1, 2*numSamples);
	if(numIntervals == 1)
		if(isArc)
			[amostra1, indexAmostra1] = sampleGaussianInterval(I, G(1), numSamples);
			if(indexAmostra1)
				branches(2*numSamples) = Ii.m;
				branchNum = 2*numSamples;
			else
				branches(numSamples + 1 : 2*numSamples) = amostra1;
				branchNum = numSamples + 1;
			end
		else
			branches(2*numSamples) = Ii.m;
			branchNum = 2*numSamples;
		end
	else
		if(isArc)
			if(numGauss == 2)
				[amostra1, indexAmostra1] = sampleGaussianInterval(I, G(1), numSamples);
				[amostra2, indexAmostra2] = sampleGaussianInterval(I, G(2), numSamples);
			else
				[amostra1, indexAmostra1] = sampleGaussianInterval(I, G(1), numSamples);
				indexAmostra2 = 1;
			end
			if(indexAmostra1 && indexAmostra2)
				branches(2*numSamples-1) = Ii(1).m;
				branches(2*numSamples) = Ii(2).m;
				branchNum = 2*numSamples-1;
			elseif(indexAmostra1 == 1 && indexAmostra2 == 0)
				branches(numSamples) = Ii(1).m;
				branches(numSamples+1 : 2*numSamples) = amostra2;
				branchNum = numSamples;
			elseif(indexAmostra1 == 0 && indexAmostra2 == 1)
				branches(numSamples : 2*numSamples-1) = amostra1;
				branches(2*numSamples) = Ii(2).m;
				branchNum = numSamples;
			else
				branches(1 : numSamples) = amostra1;
				branches(numSamples + 1 : 2*numSamples) = amostra2;
				branchNum = 1;
			end
		else
			branches(2*numSamples-1) = Ii(1).m;
			branches(2*numSamples) = Ii(2).m;
			branchNum = 2*numSamples-1;
		end
	end
end
