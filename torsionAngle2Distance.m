function d = torsionAngle2Distance(d12, d13, d23, d24, d34, Tau)

	kbar = kappai(d12, d23, d13);
	k = kappai(d24, d23, d34);
	rhobar2 = rhoi2(d12, kbar);
	rho2 = rhoi2(d24, k);
	mu = sqrt(rhobar2*rho2);
    nu = (k - kbar)*(k - kbar) + rho2 + rhobar2;
    
	d = sqrt(nu - 2*mu*cos(Tau));
end
