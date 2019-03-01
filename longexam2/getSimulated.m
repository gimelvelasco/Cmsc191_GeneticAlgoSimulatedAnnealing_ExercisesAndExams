function y_sim = getSimulated(a,x,tpl,test_func)   %Objective: The Fitness Function must be very close if not equal to 0
    if test_func == 0
        %store all a and b in x
        %x(1 to 11) is a(1 to 11)
        %x(12 to 21) is b(1 to 10)
        for k=1:numel(x)
            expr2 = a(3)*cos(a(13)*x) + a(4)*cos(a(14)*x);
            expr3 = a(7)*cos(a(17)*x) + a(8)*cos(a(18)*x) + a(9)*cos(a(19)*x);
            y_sim(k) = a(1) - a(2)*cos(a(12)*x) - expr2 - a(5)*cos(a(15)*x) - a(6)*cos(a(16)*x) - expr3 - a(10)*cos(a(20)*x) - a(11)*cos(a(21)*x);
        end
    else
        y_sim = a(1)*(x.^2) + a(2)*x + a(3);
    end
end