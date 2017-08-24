function U = legendre1_LegendreM_M_odd_fn_odd(N,m)
D = eye(N,N);

for i=0:length(D(:,1))-1
    D(i+1,i+1) = (2*i+2)*(2*i+3);
end

a = @(x) sqrt(((2*x+2)*(2*x+3)*(4*x+5))/30);

b = @(y) (m^2 - 1)*sqrt((15*(4*y+5))/(2*(2*y+2)*(2*y+3)));

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