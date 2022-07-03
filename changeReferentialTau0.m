function T0 = changeReferentialTau0(phase, T1)

	T0 = T1;
	T0.mu = phase + T0.mu;
	
	if(T0.mu < -pi)
		T0.mu = 2*pi + T0.mu;
	elseif(pi < T0.mu)
		T0.mu = -2*pi + T0.mu;
	end
end
