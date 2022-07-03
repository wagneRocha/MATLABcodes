function imprimeReslutados(idcode, X, X_ibp, X_ptaibp, I0, averageTime_ibp, stdTime_ibp, averageTime_ptaibp, stdTime_ptaibp, fiBP, fptaiBP)

	n = size(X,1);
	m = size(I0,1);
	nAA = (n -2)/5;
	
	fprintf('%s & & $%d$ & $%d$ & $%d$ & & ', idcode, nAA, n, m);
	if(fiBP)
		fprintf('$%.2f \\pm %.2f$ & ', averageTime_ibp, stdTime_ibp);
	else
		fprintf('$\\infty$ & ');
	end
	if(~isempty(X_ibp))
		satisfiesInstanceQMf(X_ibp, I0);
		[rmsd_ibp, ~] = calculaRMSD(X_ibp, X);
		[mde_ibp, lde_ibp] = calculaMDEeLDE(X_ibp, I0);
		fprintf('$%.2f$ & $%.2f$ & $%.2f$ & & ', mde_ibp, lde_ibp, rmsd_ibp);
	else
		fprintf('$--$ & $--$ & $--$ & & ');
	end
	if(fptaiBP)
		fprintf('$%.2f \\pm %.2f$ & ', averageTime_ptaibp, stdTime_ptaibp);
	else
		fprintf('$\\infty$ & ');
	end
	if(~isempty(X_ptaibp))
		satisfiesInstanceQMf(X_ptaibp, I0);
		[rmsd_ptaibp, ~] = calculaRMSD(X_ptaibp, X);
		[mde_ptaibp, lde_ptaibp] = calculaMDEeLDE(X_ptaibp, I0);
		fprintf('$%.2f$ & $%.2f$ & $%.2f$ & & ', mde_ptaibp, lde_ptaibp, rmsd_ptaibp);
%		 fprintf('$%.2f$ & $%.2f$ & $%.2f$ \\\\\n', mde_ptaibp, lde_ptaibp, rmsd_ptaibp);
	else
		fprintf('$--$ & $--$ & $--$ & & ');
%		 fprintf('$--$ & $--$ & $--$ \\\\\n');
	end

	[mde, lde] = calculaMDEeLDE(X, I0);
	fprintf('$%.2f$ & $%.2f$ & $0.00$ \\\\\n', mde, lde);
end
