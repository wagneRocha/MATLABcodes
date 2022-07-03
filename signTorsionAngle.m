function sTa = signTorsionAngle(x1, x2, x3, x4)

	v = x1 - x2;
	r = x3 - x2;
	s = x4 - x2;

	sTa = sign(dot(cross(r,v),s));
end
