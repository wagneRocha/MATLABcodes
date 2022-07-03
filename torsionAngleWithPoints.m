function tau = torsionAngleWithPoints(x1, x2, x3, x4)
	d12 = norm(x1-x2);
	d13 = norm(x1-x3);
	d14 = norm(x1-x4);
	d23 = norm(x2-x3);
	d24 = norm(x2-x4);
	d34 = norm(x3-x4);

	tau = torsionAngleWithDistances(d12, d13, d14, d23, d24, d34);
end
