function I = ajustaIntervalo(I)

	ncd = 2;

	a = I.m - I.delta;
	b = I.m + I.delta;
	if(round(a,ncd) < -pi)
		a0 = -pi;
		b0 = b;
		a1 = a;
		b1 = -pi;
		I0.m	 = (a0+b0)/2;
		I0.delta = abs(a0-b0)/2;
		I1.m	 = (a1+b1)/2;
		I1.delta = abs(a1-b1)/2;
		I1.m = 2*pi + I1.m;
		I = [I0; I1];
	elseif(round(b,ncd) > pi)
		a0 = a;
		b0 = pi;
		a1 = pi;
		b1 = b;
		I0.m	 = (a0+b0)/2;
		I0.delta = abs(a0-b0)/2;
		I1.m	 = (a1+b1)/2;
		I1.delta = abs(a1-b1)/2;
		I1.m = -2*pi + I1.m;
		I = [I0; I1];
	end
end
