function G = intersectionAmongGaussians(G1, G2, nTimesStd, myZero)

	G = [];
	for i = 1 : length(G1)
		for j = 1 : length(G2)
			res = intersectGaussians(G1(i), G2(j), nTimesStd/2, myZero);
			if(~isempty(res))
				res.delta = res.s*nTimesStd/2;
				G = [G; res];
			end
		end
	end
end
