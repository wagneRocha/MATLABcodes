function [i, exploredVertex, branches, branchNum] = backtrackingNaArvore(i, exploredVertex, branches, branchNum, numSamples)
	exploredVertex(i) = 0;
	i = i - 1;
	while(1)
		branches(i, branchNum(i)) = NaN;
		if(branchNum(i) < 2*numSamples)
			break;
		else
			exploredVertex(i) = 0;
			i = i - 1;
		end
	end
end
