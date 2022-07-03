function [mde, lde] = calculaMDEeLDE(Xsol, I)

	m = size(I,1);
	mdeVec = zeros(m,1);
	for k = 1 : m
		xi = Xsol(I(k,1),:);
		xj = Xsol(I(k,2),:);
		dijL = I(k,3);
		dijU = I(k,4);
		mdeVec(k) = max(abs(dijL - norm(xi-xj)), abs(dijU - norm(xi-xj)));
	end
	lde = max(mdeVec);
	mde = sum(mdeVec)/m;
end
