function U = legendre2_LegendreM_M_even_fn_even(N,m)
D = eye(N,N);

for i=1:length(D(:,1))
    D(i,i) = (2*i+3)*(2*i+2);
end

a = @(x) sqrt(((2*x+1)*(2*x+2)*(2*x+3)*(2*x+4)*(4*x+5))/(8*15));

b = @(y) (m^2 - 4)*sqrt((15*(4*y+5))/(2*(2*y+1)*(2*y+2)*(2*y+3)*(2*y+4)));

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
