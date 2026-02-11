

function dydt=Taylor(t,y,gamma)

    u = y(1) ;
    v = y(2) ;
    dydt = zeros(2,1);
    a_squared = (gamma-1)/2 * (1 - (u^2 + v^2));
    if a_squared < 0 
        a_squared = 1e-12;
    end
    denominator = v^2 - a_squared;
    if abs(denominator) < 1e-10
        denominator = 1e-10 * sign(denominator);
    end
    dydt(1) = v;
    dydt(2) = ((a_squared)*( (v*cot(t) + u) / ( denominator) )) - u;
    

end



