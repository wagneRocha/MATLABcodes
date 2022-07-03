function [rmsd, Xsol] = calculaRMSD(Xsol, X)

	[Xsol, X] = mesmoCentrodeMassa(Xsol, X);
	Q = melhorRotacaoNormaF(Xsol, X);
	Xsol = Xsol*Q';
	rmsd = norm(X - Xsol,'fro');
end
