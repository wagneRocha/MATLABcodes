function m = calculaM(tauFixo, tauVariavel)
	
	tau1 = tauVariavel;
	tau2 = tauVariavel - pi*sign0(tauVariavel);
	ln = min(tau1, tau2);
	lp = max(tau1, tau2);
	
	if((ln <= tauFixo) && (tauFixo <= lp))
		m = 1;
	elseif((lp <= tauFixo) && (tauFixo <= pi))
		m = 2;
	elseif((-pi <= tauFixo) && (tauFixo <= ln))
		m = 3;
	else
		printf('ERRO\n');
	end
end
