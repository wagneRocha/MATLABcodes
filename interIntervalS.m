% I = [I.mu, I.delta]
function I = interIntervalS(I1, I2, myZero)

	n = length(I1);
	m = length(I2);
	I = [];
	for i = 1 : n
		for j = 1 : m
			I = [I; interInterval(I1(i), I2(j), myZero)];
		end
	end
end
