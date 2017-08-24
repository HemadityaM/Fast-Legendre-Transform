function A = Legendre2ToChebyshev1(n)
B = zeros(n+1,n);
for k=0:length(B(:,1))-1
    for l=0:length(B(1,:))-1
        temp1 = l-k+1;
        temp2 = l+k+1;
        [Alpha1,Alpha2] = gammaFnApprox(temp1,temp2);
        if k == 0
                B(k+1,l+1) = ((2*(l+1)*(2*l+3) - 8*k^2)/(pi))*Alpha1*Alpha2*sqrt((4*l+5)/(8*(l+1)*(l+2)*(2*l+1)*(2*l+3)));
        elseif 0 < k <= (l+1)
                B(k+1,l+1) = 2*(((2*(l+1)*(2*l+3) - 8*k^2)/(pi))*Alpha1*Alpha2*sqrt((4*l+5)/(8*(l+1)*(l+2)*(2*l+1)*(2*l+3))));
                %B(k+1,l+1) = 0;
        else
            B(k+1,l+1) = 1;
        end
    end
end
A = B;
end