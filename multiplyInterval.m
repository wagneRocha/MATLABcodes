function I = multiplyInterval(I0, I1, myZero)
	
	a0aux = I0.m - I0.delta;
	b0aux = I0.m + I0.delta;
	a1aux = I1.m - I1.delta;
	b1aux = I1.m + I1.delta;

	if((a0aux < a1aux) || (abs(a0aux - a1aux) < myZero))
		b0 = b0aux;
		a1 = a1aux;
		b1 = b1aux;
	else
		b0 = b1aux;
		a1 = a0aux;
		b1 = b0aux;
	end

	if((a1 < b0) || (abs(a1 - b0) < myZero))
		a = a1;
		b = min(b0, b1);
		I.m = (a + b)/2;
		I.delta = abs(a - b)/2;
	else
		I = [];
	end
end
