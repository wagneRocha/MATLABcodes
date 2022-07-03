clc
clear all
close all

reps = 1;
myZero = 1e-6;
	
idCodes = {'1TOS', '1UAO', '1KUW', ... % 10
	   '1ID6', '1DNG', '1O53', ... % 15
	   '1DU1', '1DPK', '1HO7', ... % 20
	   '1CKZ', '1LFC', '1A11', ... % 25
	   '1HO0', '1MMC', '1D0R', ... % 30
	   '1ZWD', '1D1H', '1SPF', ... % 35
	   '1AML', '1BA4', '1C56', ... % 40
	  };

timeMax = 60;
deltaTau = 20*pi/180;
nTimesStd = 8;
numSamples = 5;

for estrutura = 1 : length(idCodes)
	idCode = idCodes{estrutura};
	fnI = ['estruturas/', idCode, '/eps_1_tau_40/', idCode, '_model1_chainA_ddgpHCorder4.dat'];
	fnX = ['estruturas/', idCode, '/eps_1_tau_40/X_', idCode, '_model1_chainA_ddgpHCorder4.dat'];
	fnphibar = ['estruturas/', idCode, '/eps_1_tau_40/phibar_', idCode, '_model1_chainA_ddgpHCorder4.dat'];
	fnpsibar = ['estruturas/', idCode, '/eps_1_tau_40/psibar_', idCode, '_model1_chainA_ddgpHCorder4.dat'];
	
	I0 = load(fnI);
	X = load(fnX);
	phibar = load(fnphibar);
	psibar = load(fnpsibar);
	m = size(I0,1);
	n = I0(m,1);
	
	AtomsOrder = DDGPvertexOrderAtoms(n);
	cliques = DDGPvertexOrder(n);
	[I, kv] = instance2MyFormat(I0, m, cliques);
	[Iu, ku] = instance2MyFormatTBP(I0, m, cliques);
	treeBranches = vertexOrderBranches(n);
	tauAngle = vertexOrderTorsionAngles(X, cliques, treeBranches, phibar, psibar);
	[tauCHA, tauNHN] = angulosFixosProteina(n, X);
	
	if(min(kv(4:n) < 3))
		fprintf('ERRO! A Instancia nao é discretizavel\n');
	else
		X0 = zeros(n, 3);
		[X0(1,:), X0(2,:), X0(3,:)] = referential_x1_x2_x3(Iu(2,1,2), Iu(3,2,2), Iu(3,1,2));
		% iTBP
		telapsed = zeros(reps,1);
		numIt_iTBP = zeros(reps,1);
		for k = 1 : reps
			tstart = tic;
			[X_iTBP, numIt_iTBP(k), numIntersec_iTBP] = iTBP(X0, Iu, n, ku, tauAngle, deltaTau, numSamples, treeBranches, myZero, nTimesStd, tauCHA, tauNHN, timeMax);
			telapsed(k) = toc(tstart);
		end
		if(not(isempty(X_iTBP)))
			iTBPflag = 1;
		else
			iTBPflag = 0;
		end
		averageTime_iTBP = mean(telapsed);
		stdTime_iTBP = std(telapsed);
		numIt_iTBP = round(mean(numIt_iTBP));
		% iBP
		telapsed = zeros(reps,1);
		numIt_iBP = zeros(reps,1);
		for k = 1 : reps
			tstart = tic;
			[X_iBP, numIt_iBP(k), numIntersec_iBP] = iBP(X0, I, n, kv, tauAngle, deltaTau, numSamples, treeBranches, myZero, nTimesStd, timeMax);
			telapsed(k) = toc(tstart);
		end
		if(not(isempty(X_iBP)))
			iBPflag = 1;
		else
			iBPflag = 0;
		end
		averageTime_iBP = mean(telapsed);
		stdTime_iBP = std(telapsed);
		numIt_iBP = round(mean(numIt_iBP));
		
		imprimeReslutados(idCode, X, X_iBP, X_iTBP, I0, averageTime_iBP, stdTime_iBP, averageTime_iTBP, stdTime_iTBP, iBPflag, iTBPflag)
	end
end
