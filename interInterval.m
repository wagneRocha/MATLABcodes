% I.mu \in [-pi,pi]
function I = interInterval(I1, I2, myZero)
	I1 = ajustaIntervalo(I1);
	I2 = ajustaIntervalo(I2);
	I = [];
	for i = 1 : length(I1)
		for j = 1 : length(I2)
			I = [I; multiplyInterval(I1(i), I2(j), myZero)];
		end
	end
end
