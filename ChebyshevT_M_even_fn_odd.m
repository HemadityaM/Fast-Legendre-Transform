function interp = ChebyshevT_M_even_fn_odd(coeffs,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N-1
        interp = interp + coeffs(i+1)*chebyshevT(2*i + 1,x);
    end
end