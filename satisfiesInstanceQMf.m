function satisfiesInstanceQMf(X, I)
	
	[numberOfMeasures,~] = size(I);
	flags = ones(numberOfMeasures,1);
	
	l = 0;
	for k = 1 : numberOfMeasures
		xi = X(I(k,1),:);
		xj = X(I(k,2),:);
		
		d  = round(norm(xi - xj),4);
		dl = round(I(k,3),4);
		du = round(I(k,4),4);

		if((dl <= d) && (d <= du))
			flags(k) = 0;
		else
			l = l + 1;
		end
	end
	if(sum(flags) == 0)
%		 fprintf('teste: X Satisfies the Instance\n')
	else
		fprintf('teste: X does not Satisfies the Instance: %d wrong edges were detected\n', l)
		vec = find(flags == 1)';
		for k = 1 : length(vec)
			[I(vec(k),:) 0 norm(X(I(vec(k),1),:) - X(I(vec(k),2),:))]
		end
	end
end
