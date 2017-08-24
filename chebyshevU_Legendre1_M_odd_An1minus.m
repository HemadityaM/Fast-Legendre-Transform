function B = chebyshevU_Legendre1_M_odd_An1minus(n)
A = zeros(n,n);
for l=0:length(A(:,1))-1
    for k=0:length(A(1,:))-1
        temp1 = k-l;
        temp2 = (2*k+2*l+3)/2;
        [Alpha1,Alpha2] = gammaFnApprox(temp1,temp2);
        if ge(k,l)
                A(l+1,k+1) = -((4*l+5)/(2*(2*k+2*l+5)*(2*k-2*l-1)))*Alpha1*Alpha2*sqrt((4*(l+1)*(2*l+3))/(4*l+5));
                
        else
            A(l+1,k+1) = 0;
        end
    end
end
B = A;
end