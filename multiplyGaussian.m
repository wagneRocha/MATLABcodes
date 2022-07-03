% g :   g.mu
%   :   g.s
%   :   g.A
function g = multiplyGaussian(g1, g2)
	
	g.mu = (g1.mu*g2.s*g2.s + g2.mu*g1.s*g1.s)/(g1.s*g1.s + g2.s*g2.s);
	g.s = sqrt((g1.s*g1.s*g2.s*g2.s)/(g1.s*g1.s + g2.s*g2.s));
	g.A = g1.A*g2.A*exp(-0.5*(g1.mu - g2.mu)*(g1.mu - g2.mu)/(g1.s*g1.s + g2.s*g2.s));
end
