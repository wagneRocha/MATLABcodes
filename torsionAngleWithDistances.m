function tau = torsionAngleWithDistances(d12, d13, d14, d23, d24, d34)

	kbar = kappai(d12, d23, d13);
	k = kappai(d24, d23, d34);
	rhobar2 = rhoi2(d12, kbar);
	rho2 = rhoi2(d24, k);
	
	nu = (k - kbar)*(k - kbar) + rhobar2 + rho2;
	mu = sqrt(rhobar2*rho2);
    
	dmin = sqrt(nu-2*mu);
	dmax = sqrt(nu+2*mu);
    
	if((dmin < d14) && (d14 < dmax))
		cost = round((nu - d14*d14)/(2*mu), 12);
	elseif(d14 <= dmin)
		cost = 1;
	elseif(dmax <= d14)
		cost = -1;
	end
	
	tau = acos(cost);
end
