function nAA = numeroDeAminoAcidos(n)
   
	caso = mod(n,5);
	if(caso == 0)
		% nao tem H3 e H2
		desloca = 0;
	elseif(caso == 1)
		% tem H2 ou H3
		desloca = 1;
	elseif(caso == 2)
		% tem H2 e H3
		desloca = 2;
	else
		fprintf('ERRO\n');
	end
	nAA = (n-desloca)/5;
end
