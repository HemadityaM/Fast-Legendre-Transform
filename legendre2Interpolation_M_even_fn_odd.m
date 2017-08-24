function interp = legendre2Interpolation_M_even_fn_odd(coeffs,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N-1
        P = legendre(2*i + 3, x,'norm');
        interp = interp + coeffs(i+1)*P(3,:);
    end
end