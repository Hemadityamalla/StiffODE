function f = temperature(t)
%Write the tanh function for temperature here
    if t < 4e-5
        f = 1200;
    else
        f = 2700*tanh(0.5e5*(t - 3.95e-5));
    end
end
