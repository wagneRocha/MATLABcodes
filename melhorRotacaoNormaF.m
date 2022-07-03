function Q = melhorRotacaoNormaF(Y,X)

	[U, ~, V] = svd(X'*Y);
	Q = U*V';
end
