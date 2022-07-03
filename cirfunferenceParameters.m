function [a,b,c] = cirfunferenceParameters(x1, x2, x3, d24, d34)
	
	d12 = norm(x1-x2);
	d13 = norm(x1-x3);
	d23 = norm(x2-x3);

	rhat = (x3 - x2)/norm(x3 - x2);
	v = x1 - x2;
	kbar = kappai(d12, d23, d13);
	k = kappai(d24, d23, d34);
	rhobar2 = rhoi2(d12, kbar);
	rho2 = rhoi2(d24, k);
	mu = sqrt(rhobar2*rho2);
	
	a = x2 + k*rhat;
	b = (rho2/mu)*(v - kbar*rhat);
	c = (rho2/mu)*cross(rhat,v);
end
