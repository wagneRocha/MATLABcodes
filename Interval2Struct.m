function T = Interval2Struct(I, myZero)
	
	n = size(I,1);
	for k = 1 : n
		T(k).m = (I(k,1) + I(k,2))/2;
		T(k).delta = abs(I(k,1) - I(k,2))/2;
	end
	
	T11 = T(1).m - T(1).delta;
	Tend2 = T(end).m + T(end).delta;
	if(((T11 < -pi) || (abs(T11 + pi) < myZero)) && ((pi < Tend2) || (abs(Tend2 - pi) < myZero)))
		a = T(end).m - T(end).delta - 2*pi;
		b = T(1).m + T(1).delta;
		T(1).m = (a + b)/2;
		T(1).delta = abs(a - b)/2;
		T(end) = [];
	end
end
