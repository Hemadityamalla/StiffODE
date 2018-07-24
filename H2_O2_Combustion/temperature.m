function f = temperature(t)
%Write the tanh function for temperature here
    if t < 3.95e-5
        f = 273;
    else
        f = 2700*tanh(5e4*(t - 3.95e-5));
    end
end
