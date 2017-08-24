function U = legendre2_LegendreM_M_even_fn_odd(N,m)

D = eye(N,N);

for i=0:length(D(:,1))-1
    D(i+1,i+1) = (2*i+3)*(2*i+4);
end

a = @(x) sqrt(((2*x+2)*(2*x+3)*(2*x+4)*(2*x+5)*(4*x+7))/(14*15));

b = @(y) (m^2 - 4)*sqrt((7*15*(4*y+7))/(8*(2*y+2)*(2*y+3)*(2*y+4)*(2*y+5)));

S = zeros(N,N);

for j=0:length(S(1,:))-1
    for k=0:length(S(:,1))-1
        if ge(k,j) 
            S(j+1,k+1) = a(j)*b(k);
        else
            S(j+1,k+1) = a(k)*b(j);
            %S(j+1,k+1) = S(j+1,k+1);
        end      
    end
end
[U,V] = eig(D+S);

end
