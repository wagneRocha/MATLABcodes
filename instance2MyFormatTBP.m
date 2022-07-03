function [I, kv] = instance2MyFormatTBP(I0, m, vertexOrder)

	[Iaux, kv] = instance2MyFormat(I0, m, vertexOrder);
	n = length(kv);
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
	nmeU = max(kv);
	
	I = zeros(n,nmeU + 3,3);
	% I = [atomo_i, [atomo_j (predecessor adjacente), dij, Dij]]
	
	I(1,1:3,:) = Iaux(1,1:3,:);
	I(2,1:3,:) = Iaux(2,1:3,:);
	I(3,1:3,:) = Iaux(3,1:3,:);
	I(4,1:3,:) = Iaux(4,1:3,:);
	I(5,1:3,:) = Iaux(5,1:3,:);
	if(caso == 0)
		I(5,4:nmeU,:) = Iaux(5,4:nmeU,:);
	elseif(caso == 1)
		I(6,1:3,:) = Iaux(6,1:3,:);
		I(5,4:nmeU,:) = Iaux(5,4:nmeU,:);
	else
		I(6,1:3,:) = Iaux(6,1:3,:);
		I(7,1:3,:) = Iaux(7,1:3,:);
		I(6,4:nmeU,:) = Iaux(6,4:nmeU,:);
	end
	
	for i = 2 : nAA
		ai1 = desloca + 5*(i-1) + 1;
		I(ai1+0,1:3,:) = Iaux(ai1+0,1:3,:);
		I(ai1+1,1:3,:) = Iaux(ai1+1,1:3,:);
		I(ai1+2,1:3,:) = Iaux(ai1+2,1:3,:);
		I(ai1+3,1:3,:) = Iaux(ai1+3,1:3,:);
		I(ai1+4,1:3,:) = Iaux(ai1+4,1:3,:);
		
		I(ai1+0,4:nmeU+3,:) = Iaux(ai1+1,:,:);
		I(ai1+3,4:nmeU+3,:) = Iaux(ai1+4,:,:);
	end
	
	for v = 1 : n
		kv(v,1) = length(find(I(v,:,1) > 0));
	end
end
