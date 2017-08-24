function interp = legendreMInterpolation_M_odd_fn_even(coeffs,m,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N - (m-1)/2
        P = legendre(2*i + m ,x,'norm');
        interp = interp + coeffs(i+1)*P(m+1,:);
    end
end