function T = juntaIntervalos(I, myZero)
	n = length(I);
	Taux = zeros(n,2);
	for k = 1 : n
		Taux(k,:) = [I(k).m - I(k).delta, I(k).m + I(k).delta];
	end
	[Taux,~,~] = unique(Taux, 'rows');
	T = Taux(1,:);
	l = 1;
	n = size(Taux,1);
	for k = 1 : n-1
		if(abs(Taux(k,2) - Taux(k+1,1)) < myZero)
			T(l,2) = Taux(k+1,2);
		else
			T = [T; Taux(k+1,:)];
			l = l + 1; 
		end
	end
end
