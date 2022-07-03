function vertexOrder = DDGPvertexOrder(n)

	caso = mod(n,5);
	if(caso == 0)
		% nao tem H3 e H2
		desloca = 2;
	elseif(caso == 1)
		% tem H2 ou H3
		desloca = 1;
	elseif(caso == 2)
		% tem H2 e H3
		desloca = 0;
	else
		fprintf('ERRO\n');
	end

	vertexOrder = zeros(n,4);

	vertexOrderRes1 = zeros(7-desloca,4);
	if(caso == 0)
		for i = 1 : 4
			vertexOrderRes1(i,:) = [max(i,0), max(i-1,0), max(i-2,0), max(i-3,0)];
		end
		vertexOrderRes1(5,:) = [5, 1, 2, 4];
	else
		for i = 1 : 7-desloca
			vertexOrderRes1(i,:) = [max(i,0), max(i-1,0), max(i-2,0), max(i-3,0)];
		end
	end
	% DDGP hc Order:{H3, H2, H1, N1, CA1, HA1, C1, // N1, CA1, C1, HA1, H1,
	%		...,
	%		Ni, HZi, CAi, Ci, HAi
	%		...,
	%		Np, HZp, CAp, Cp, HAp}
	% -------------------- 2nd to p-th Residue --------------------
	if(caso == 0)
		vertexOrderRes2 = [ 6, 3, 2, 1
				    7, 6, 3, 2
				    8, 7, 6, 3
				    9, 8, 6, 3
				    10, 9, 8, 6];
	else
		vertexOrderRes2 = [ 8, 7, 5, 4
				    9, 8, 7, 5
				    10, 9, 8, 7
				    11, 10, 8, 7
				    12, 11, 10, 8]-desloca;
	end

	vertexOrder(1:12-desloca,:) = [vertexOrderRes1; vertexOrderRes2];

	for i = 13-desloca : n
		switch mod(i-(12-desloca),5)
			case 1 % N
				vertexOrder(i,:) = [i, i-2, i-3, i-5];
			case 2 % HZ
				vertexOrder(i,:) = [i, i-1, i-3, i-4];
			case 3 % CA
				vertexOrder(i,:) = [i, i-1, i-2, i-4];
			case 4 % C
				vertexOrder(i,:) = [i, i-1, i-3, i-5];
			case 0 % HA
				vertexOrder(i,:) = [i, i-1, i-2, i-4];
		end
	end
end
