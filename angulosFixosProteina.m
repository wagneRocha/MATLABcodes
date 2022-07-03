function [tauCHA, tauNHN] = angulosFixosProteina(n, X)

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
	nAA = numeroDeAminoAcidos(n);
	
	tauCHA = zeros(nAA,1);
	if(caso == 0)
		% nao tem H3 e H2
		x1 = X(3,:);
		x2 = X(1,:);
		x3 = X(2,:);
		x4 = X(4,:);
	elseif(caso == 1)
		% tem H2 ou H3
		x1 = X(6,:);
		x2 = X(3,:);
		x3 = X(4,:);
		x4 = X(5,:);
	else
		% tem H2 e H3
		x1 = X(7,:);
		x2 = X(4,:);
		x3 = X(5,:);
		x4 = X(6,:);
	end
	tauCHA(1) = torsionAngleWithPoints(x1, x2, x3, x4)*signTorsionAngle(x1, x2, x3, x4);
	for k = 2 : nAA
		ai1 = desloca + 5*(k-1)+1;
		x1 = X(ai1+3,:);
		x2 = X(ai1+0,:);
		x3 = X(ai1+2,:);
		x4 = X(ai1+4,:);
		tauCHA(k) = torsionAngleWithPoints(x1, x2, x3, x4)*signTorsionAngle(x1, x2, x3, x4);
	end
	
	tauNHN = zeros(nAA,1);
	if(caso == 0)
		% nao tem H3 e H2
		x1 = X(6,:);
		x2 = X(2,:);
		x3 = X(3,:);
		x4 = X(7,:);
	elseif(caso == 1)
		% tem H2 ou H3
		x1 = X(7,:);
		x2 = X(4,:);
		x3 = X(6,:);
		x4 = X(8,:);
	else
		% tem H2 e H3
		x1 = X(8,:);
		x2 = X(5,:);
		x3 = X(7,:);
		x4 = X(9,:);
	end
	tauNHN(2) = torsionAngleWithPoints(x1, x2, x3, x4)*signTorsionAngle(x1, x2, x3, x4);
	for k = 3 : nAA
		ai1 = desloca + 5*(k-1)+1;
		x1 = X(ai1+0,:);
		x2 = X(ai1-3,:);
		x3 = X(ai1-2,:);
		x4 = X(ai1+1,:);
		tauNHN(k) = torsionAngleWithPoints(x1, x2, x3, x4)*signTorsionAngle(x1, x2, x3, x4);
	end
end
