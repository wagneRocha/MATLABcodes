function AtomsOrder = DDGPvertexOrderAtoms(n)
	caso = mod(n,5);
	if(caso == 0)
		% nao tem H3 e H2
		AtomsOrder = {'N', 'CA', 'C', 'HA', 'H'};
		desloca = 2;
	elseif(caso == 1)
		% tem H2 ou H3
		AtomsOrder = {'H2', 'H', 'N', 'CA', 'HA', 'C'};
		desloca = 1;
	elseif(caso == 2)
		% tem H2 e H3
		AtomsOrder = {'H3', 'H2', 'H', 'N', 'CA', 'HA', 'C'};
		desloca = 0;
	else
		fprintf('ERRO\n');
	end

	% DDGP hc Order: {H3, H2, H1, N1, CA1, HA1, C1,
	%				 ...,
	%				 Ni, HZi, CAi, Ci, HAi
	%				 ...,
	%				 Np, HZp, CAp, Cp, HAp}
	% -------------------- 2nd to p-th Residue --------------------
	for i = 1 : (n-7-desloca)/5
		AtomsOrder = [AtomsOrder, {'N', 'H', 'CA', 'C', 'HA'}];
	end
	AtomsOrder = AtomsOrder';
end
