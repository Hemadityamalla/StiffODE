function f = dTdt(t)
    if t < 3.95e-5
        f = 0;
    else
        f = 1.35e8*sech(1.975 - 5e4*t)^2;
    end
end