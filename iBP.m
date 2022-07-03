function [X, numIt, numIntersec] = iBP(X, I, n, kv, tauAngle, deltaTau, numSamples, isArc, myZero, nTimesStd, timeMax)
	
	A = zeros(n,3);
	B = A;
	C = A;
	numIt = 0;
	numIntersec = 0;
	
	branches = NaN(n,2*numSamples);
	exploredVertex = zeros(n);
	branchNum = ones(n,1);
	
	i = 4;
	iniTime = clock;
	while((i <= n) && (etime(clock, iniTime) < timeMax))
%		 fprintf('iBP: layer %d\n', i);
		if(exploredVertex(i) == 0)
			% jn < jn-1 < ... < j3 < j2 < j1
			Jn = zeros(3,1);
			for j = 1 : kv(i)
				Jn(j) = I(i,j,1);
			end
			u1 = X(Jn(1),:);
			u2 = X(Jn(2),:);
			u3 = X(Jn(3),:);
			diu1 = I(i, 1, 2);
			diu2 = I(i, 2, 2);
			diu3L = I(i, 3, 2);
			diu3U = I(i, 3, 3);
%			 fprintf('terceira distancia\n');
			tau0l = torsionAngleWithPointsAndDistances(u3, u2, u1, diu3L, diu2, diu1);
			tau0u = torsionAngleWithPointsAndDistances(u3, u2, u1, diu3U, diu2, diu1);
			
			clear G0 T0
			G0(1).mu	= (tau0u + tau0l)/2;
			G0(1).s	 = (tau0u - tau0l)/nTimesStd;
			G0(1).delta = G0(1).s*nTimesStd/2;
			G0(1).A	 = 1;
			
			G0(2) = G0(1);
			G0(2).mu = -G0(2).mu;
			
			T0(1).m	 = G0(1).mu;
			T0(1).delta = G0(1).delta;
			T0(2).m	 = G0(2).mu;
			T0(2).delta = G0(2).delta;
			
			G = [];
			T = [];
			if(~isnan(tauAngle(i)))
%				 fprintf('angulo de torcao a priori\n');
				G.mu = tauAngle(i);
				G.A  = 1;
				if(isArc(i))
					G.s	 = 2*deltaTau/nTimesStd;
					G.delta = deltaTau;
				else
					G.s	 = 0;
					G.delta = 0;
				end
				T.m	 = G.mu;
				T.delta = G.delta;
			end

			if(isempty(G))
				G = G0;
				T = T0;
			end
			
			% jn < jn-1 < ... < j4
			Hn = zeros(kv(i)-3,1);
			for j = 1 : kv(i)-3
				Hn(j) = I(i,j+3,1);
			end
			for k = 4 : kv(i)
%				 fprintf('distancia extra %d\n', k);
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
				
				T1pk.mu	= (tauKu + tauKl)/2;
				T1pk.s	 = abs(tauKu - tauKl)/nTimesStd;
				T1pk.delta = T1pk.s*nTimesStd/2;
				T1pk.A	 = 1;
				
				T1nk = T1pk;
				T1nk.mu = -T1nk.mu;
				
				clear Gk
				Gk(1) = changeReferentialTau0(phasek, T1pk);
				Gk(2) = changeReferentialTau0(phasek, T1nk);
				
				clear Tk
				Tk(1).m	 = Gk(1).mu;
				Tk(1).delta = Gk(1).delta;
				Tk(2).m	 = Gk(2).mu;
				Tk(2).delta = Gk(2).delta;
				
				G = intersectionAmongGaussians(G, Gk, nTimesStd, myZero);
				T = interIntervalS(T, Tk, myZero);
				numIntersec = numIntersec + 1;
			end
			
			if((~isempty(G)) && (~isempty(T)))
				[branches(i,:), branchNum(i)] = sampleInterval(T, G, isArc(i), numSamples, myZero);
				exploredVertex(i) = 1;
				[A(i,:), B(i,:), C(i,:)] = cirfunferenceParameters(u3, u2, u1, diu2, diu1);
			else
				% backtracking na arvore
				[i, exploredVertex, branches, branchNum] = backtrackingNaArvore(i, exploredVertex, branches, branchNum, numSamples);
				if(i == 3)
%					 fprintf('ibp: ERRO! Nao encontrei a solucao da instancia\n');
					break;
				end
				continue;
			end
		else
			branches(i, branchNum(i)) = NaN;
			branchNum(i) = branchNum(i) + 1;
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
%				 fprintf('ibp: ERRO! Nao encontrei a solucao da instancia\n');
				break;
			end
			continue;
		end
	end
	if(i <= n)
		X = [];
	end
end
