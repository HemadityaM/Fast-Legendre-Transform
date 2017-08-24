function B = chebyshevT_Legendre2_M_even_An2plus(n)
A = zeros(n,n+1);
for l=0:length(A(:,1))-1
    for k=0:length(A(1,:))-1
        temp1 = k-l-1;
        temp2 = (2*k+2*l+1)/2;
        [Alpha1,Alpha2] = gammaFnApprox(temp1,temp2);
        if ge(k,l+1)
                A(l+1,k+1) = (4 + ((2*k*(6*(l+1)*(2*l+3) - 2*(2*k-1)*(2*k+1)))/((2*k+2*l+3)*(2*k-2*l-3)))*Alpha1*Alpha2)*...
                sqrt((4*l+5)/(8*(l+1)*(l+2)*(2*l+1)*(2*l+3)));
                
        else
            A(l+1,k+1) = 4*sqrt((4*l+5)/(8*(l+1)*(l+2)*(2*l+1)*(2*l+3)));
        end
    end
end
B = A;
end