function [x1,x2,x3] = referential_x1_x2_x3(d12, d13, d23)
	
	cosTh = (d12*d12 + d23*d23 - d13*d13)/(2*d12*d23);
	sinTh = sqrt(1 - cosTh*cosTh);
	
	x1 = [-d12*sinTh, d12*cosTh, 0];
	x2 = [0, 0, 0];
	x3 = [0, d23, 0];
end
