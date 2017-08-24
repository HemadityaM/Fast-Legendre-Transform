function interp = legendreMInterpolation_M_even_fn_odd(coeffs,m,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N-1 %- m/2 - 1
        P = legendre(2*i + m + 1,x,'norm');
        interp = interp + coeffs(i+1)*P(m+1,:);
    end
end
