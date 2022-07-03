function G = intersectGaussians(g1, g2, nTimesStd, myZero)

	if(abs(g1.mu - g2.mu) > pi)
		if(sign(g1.mu) < 0)
			g1.mu = 2*pi + g1.mu;
		else
			g2.mu = 2*pi + g2.mu;
		end
	end
	
	if(g1.mu - g2.mu > 1e-4) % caso 1: g2.mu < g1.mu
		
		if((g1.mu - nTimesStd*g1.s) < (g2.mu + nTimesStd*g2.s))
			G = multiplyGaussian(g1, g2);
			if(G.mu > pi)
				G.mu = -2*pi + G.mu;
			end
		else
			G = [];
		end
	elseif(g2.mu - g1.mu > 1e-4) % caso 2: g1.mu < g2.mu
		
		if((g2.mu - nTimesStd*g2.s) < (g1.mu + nTimesStd*g1.s))
			G = multiplyGaussian(g1, g2);
			if(G.mu > pi)
				G.mu = -2*pi + G.mu;
			end
		else
			G = [];
		end
	else % caso 3: g1.mu = g2.mu
		G = g1;
	end
end
