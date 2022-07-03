function [Xbar, Ybar] = mesmoCentrodeMassa(X,Y)
	
	[n,~] = size(X);
	MX = ones(n,3);
	MY = MX;
	CMX = sum(X)/n;
	CMY = sum(Y)/n;
	MX(:,1) = CMX(1)*MX(:,1);
	MX(:,2) = CMX(2)*MX(:,2);
	MX(:,3) = CMX(3)*MX(:,3);
	MY(:,1) = CMY(1)*MY(:,1);
	MY(:,2) = CMY(2)*MY(:,2);
	MY(:,3) = CMY(3)*MY(:,3);
	Xbar = X - MX;
	Ybar = Y - MY;
end
