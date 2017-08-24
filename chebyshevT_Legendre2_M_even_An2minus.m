function B = chebyshevT_Legendre2_M_even_An2minus(n)
A = zeros(n-1,n);
for l=0:length(A(:,1))-1
    for k=0:length(A(1,:))-1
        temp1 = k-l-1;
        temp2 = (2*k+2*l+3)/2;
        [Alpha1,Alpha2] = gammaFnApprox(temp1,temp2);
        if ge(k,l+1)
                A(l+1,k+1) = (4 + (((2*k+1)*(6*(l+2)*(2*l+3) - 8*k*(k+1)))/((2*k+2*l+5)*(2*k-2*l-3)))*Alpha1*Alpha2)*...
                sqrt((4*l+7)/(8*(l+1)*(l+2)*(2*l+3)*(2*l+5)));
                
        else
            A(l+1,k+1) = 4*sqrt((4*l+7)/(8*(l+1)*(l+2)*(2*l+3)*(2*l+5)));
        end
    end
end
B = A;
end