function nextT = cooling_implementation(ci,ca,cb,cooling_ratio,sched,T)
	if ci == 1	%Geometric Cooling. it is good at looking for the coefficients but it does have a quite negative effect since it easily converges to a particular solution
		a = ca;
		nextT = a*T;
	elseif ci == 2	%The Another Cooling method w/ a small b. this one is never used in the program since it takes up way too much time to cool down specially when T is below 1.
		b = cb;
		nextT = T/(1 + (b*T));
	else 	 	%Exponential Cooling. this cooling implementation is fast and is much more accurate when the cooling ratio is closer to 1
		nextT = T*(cooling_ratio)^sched;
	end
end