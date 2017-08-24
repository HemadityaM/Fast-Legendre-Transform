function interp = ChebyshevU_M_odd_fn_odd(coeffs,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N - 1
        interp = interp + coeffs(i+1)*chebyshevU(2*i + 1,x).*sqrt(1 - x.^2);
    end
end