function c = calculaC(xi, m, tauVariavel)
	
	if(abs(xi) < pi)
		c = 0;
	elseif(abs(xi) > pi)
		c = 2*pi*((-1)^(m+1))*((sign0(tauVariavel))^(mod(m,3)));
	elseif(abs(xi) == pi)
		c = pi - xi;
	else
		printf('ERRO\n');
	end
end
