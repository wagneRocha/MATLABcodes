function tau = torsionAngleWithPointsAndDistances(x1, x2, x3, d14, d24, d34)

	d12 = norm(x1-x2);
	d13 = norm(x1-x3);
	d23 = norm(x2-x3);
    
	tau = torsionAngleWithDistances(d12, d13, d14, d23, d24, d34);
end
