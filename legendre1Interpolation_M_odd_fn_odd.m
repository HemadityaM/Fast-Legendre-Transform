function interp = legendre1Interpolation_M_odd_fn_odd(coeffs,x)
    N = length(coeffs);
    interp = 0.0;
    for i=0:N-1
        P = legendre(2*i + 2, x,'norm');
        interp = interp + coeffs(i+1)*P(2,:);
    end
end