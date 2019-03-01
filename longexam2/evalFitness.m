function nrmsd = evalFitness(y_sim,y,xmax,xmin)   %Objective: The Fitness Function must be very close if not equal to 0
    oexpr = 0;
    iexpr = ((y_sim - y).^2)/(xmax-xmin);
    for j=1:numel(y)
        oexpr = oexpr + iexpr(j);
    end
    nrmsd = sqrt(oexpr);

end