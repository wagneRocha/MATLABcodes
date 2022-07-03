function [X, numIt, numIntersec] = iTBP(X, I, n, kv, tauAngle, deltaTau, numSamples, isArc, myZero, nTimesStd, tauCHA, tauNHN, timeMax)
	
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
	
	A = zeros(n,3);
	B = A;
	C = A;
	numIt = 0;
	numIntersec = 0;
	
	branches = nan(n,2*numSamples);
	exploredVertex = zeros(n);
	branchNum = ones(n,1);
	
	i = 4;
	iniTime = clock;
	while((i <= n) && (etime(clock, iniTime) < timeMax))
%		 fprintf('ptaBP: layer %d\n', i);
		% j3 < j2 < j1
		Jn = zeros(3,1);
		for j = 1 : 3
			Jn(j) = I(i,j,1);
		end
		u1 = X(Jn(1),:);
		u2 = X(Jn(2),:);
		u3 = X(Jn(3),:);
		diu1 = I(i, 1, 2);
		diu2 = I(i, 2, 2);
		if(~isArc(i))
			[A(i,:), B(i,:), C(i,:)] = cirfunferenceParameters(u3, u2, u1, diu2, diu1);
			branchNum(i) = 2*numSamples;
			branches(i, 2*numSamples) = tauAngle(i);
			tau = branches(i, branchNum(i));
			X(i,:) = A(i,:) + B(i,:)*cos(tau) + C(i,:)*sin(tau);
			numIt = numIt + 1;
			i = i + 1;
			continue;
		else
			if(exploredVertex(i) == 0)
				if(i == 6)
					diu3L = I(i, 3, 2);
					diu3U = I(i, 3, 3);
					
					tau0l = torsionAngleWithPointsAndDistances(u3, u2, u1, diu3L, diu2, diu1);
					tau0u = torsionAngleWithPointsAndDistances(u3, u2, u1, diu3U, diu2, diu1);
					
					clear G T
					G(1).mu	= (tau0u + tau0l)/2;
					G(1).s	 = (tau0u - tau0l)/nTimesStd;
					G(1).delta = G(1).s*nTimesStd/2;
					G(1).A	 = 1;
					
					G(2) = G(1);
					G(2).mu = -G(2).mu;
					
					T(1).m	 = G(1).mu;
					T(1).delta = G(1).delta;
					T(2).m	 = G(2).mu;
					T(2).delta = G(2).delta;
					
					% jn < jn-1 < ... < j4 para Hk
					Hn = zeros(kv(i)-3,1);
					for j = 1 : kv(i)-3
						Hn(j) = I(i,j+3,1);
					end
					for k = 4 : kv(i)
						uk = X(Hn(k-3),:);
						dukiL = I(i, k, 2);
						dukiU = I(i, k, 3);
						
						phasek = torsionAngleWithPoints(u3, u2, u1, uk)*signTorsionAngle(u3, u2, u1, uk);
						
						tauKl = torsionAngleWithPointsAndDistances(uk, u2, u1, dukiL, diu2, diu1);
						tauKu = torsionAngleWithPointsAndDistances(uk, u2, u1, dukiU, diu2, diu1);
						
						if(((tauKl == pi) && (tauKu == pi)) || ((tauKl == 0) && (tauKu == 0)))
							T = [];
							G = [];
							break;
						end
				
						tHkHz.mu	= (tauKu + tauKl)/2;
						tHkHz.s	 = abs(tauKu - tauKl)/nTimesStd;
						tHkHz.delta = tHkHz.s*nTimesStd/2;
						tHkHz.A	 = 1;
						
						tauKp = tHkHz;
						tauKn = tHkHz;
						tauKn.mu = -tauKn.mu;
						
						clear Gk
						Gk(1) = changeReferentialTau0(phasek, tauKp);
						Gk(2) = changeReferentialTau0(phasek, tauKn);
						
						clear Tk
						Tk(1).m	 = Gk(1).mu;
						Tk(1).delta = Gk(1).delta;
						Tk(2).m	 = Gk(2).mu;
						Tk(2).delta = Gk(2).delta;
						
						G = intersectionAmongGaussians(G, Gk, nTimesStd, myZero);
						T = interIntervalS(T, Tk, myZero);
						numIntersec = numIntersec + 1;
					end
				else
					clear G T
					G.mu	= tauAngle(i);
					G.s	 = 2*deltaTau/nTimesStd;
					G.delta = deltaTau;
					G.A	 = 1;
					T.m	 = G.mu;
					T.delta = G.delta;
					
					% jn < jn-1 < ... < j4 < j3 < j2 < j1 para Hk
					Hn = zeros(max(kv)-3,1);
					for j = 1 : kv(i)-3
						Hn(j) = I(i, j+3, 1);
					end
					du1Hz = I(i, 5, 2);
					du2Hz = I(i, 6, 2);
					
					resAtom_i = ceil((i-desloca)/5);
					if(mod(i-desloca,5) == 1)
						tauVariavel = tauNHN(resAtom_i);
					elseif(mod(i-desloca,5) == 4)
						tauVariavel = tauCHA(resAtom_i);
					end
					
					for k = 1 : kv(i)-6
						Hk = X(Hn(k+3),:);
						dHkHzL = I(i, k+6, 2);
						dHkHzU = I(i, k+6, 3);
						
						tauFixo = torsionAngleWithPoints(u3, u2, u1, Hk)*signTorsionAngle(u3, u2, u1, Hk);
						m = calculaM(tauFixo, tauVariavel);
						tau_ijk = tauFixo - tauVariavel;
						
						tauHkHzl = torsionAngleWithDistances(norm(Hk-u2), norm(Hk-u1), dHkHzL, norm(u1-u2), du2Hz, du1Hz);
						tauHkHzu = torsionAngleWithDistances(norm(Hk-u2), norm(Hk-u1), dHkHzU, norm(u1-u2), du2Hz, du1Hz);
						
						if(((tauHkHzl == pi) && (tauHkHzu == pi)) || ((tauHkHzl == 0) && (tauHkHzu == 0)))
							T = [];
							G = [];
							break;
						end
						
						tHkHz.mu	= (tauHkHzl + tauHkHzu)/2;
						tHkHz.s	 = abs(tauHkHzl - tauHkHzu)/nTimesStd;
						tHkHz.delta = tHkHz.s*nTimesStd/2;
						tHkHz.A	 = 1;
						
						clear Gk
						Gk(1) = tHkHz;
						Gk(2) = tHkHz;
						
						xip = tau_ijk + tHkHz.mu;
						xin = tau_ijk - tHkHz.mu;
						cp = calculaC(xip, m, tauVariavel);
						cn = calculaC(xin, m, tauVariavel);
						
						Gk(1).mu = cp + xip;
						Gk(2).mu = cn + xin;
						
						clear Tk
						Tk(1).m	 = Gk(1).mu;
						Tk(1).delta = Gk(1).delta;
						Tk(2).m	 = Gk(2).mu;
						Tk(2).delta = Gk(2).delta;
						
						G = intersectionAmongGaussians(G, Gk, nTimesStd, myZero);
						T = interIntervalS(T, Tk, myZero);
						numIntersec = numIntersec + 1;
					end
				end
				
				if((~isempty(G)) && (~isempty(T)))
					[branches(i,:), branchNum(i)] = sampleInterval(T, G, isArc(i), numSamples, myZero);
					exploredVertex(i) = 1;
					[A(i,:), B(i,:), C(i,:)] = cirfunferenceParameters(u3, u2, u1, diu2, diu1);
				else
					% backtracking na arvore
					[i, exploredVertex, branches, branchNum] = backtrackingNaArvore(i, exploredVertex, branches, branchNum, numSamples);
					if(i == 3)
%						 fprintf('ptaBP: ERRO! Nao encontrei a solucao da instancia\n');
						break;
					end
					continue;
				end
			else
				branches(i, branchNum(i)) = NaN;
				branchNum(i) = branchNum(i) + 1;
			end
		end
		
		if(branchNum(i) <= 2*numSamples) % (o ultimo ramo nao foi explorado)
			tau = branches(i, branchNum(i));
			X(i,:) = A(i,:) + B(i,:)*cos(tau) + C(i,:)*sin(tau);
			numIt = numIt + 1;
			i = i + 1;
			continue;
		else
			% backtracking na arvore
			[i, exploredVertex, branches, branchNum] = backtrackingNaArvore(i, exploredVertex, branches, branchNum, numSamples);
			if(i == 3)
%				 fprintf('ptaBP: ERRO! Nao encontrei a solucao da instancia\n');
				break;
			end
			continue;
		end
	end
	if(i <= n)
		X = [];
	end
end
