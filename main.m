clear all;
clc;

N = 128;

m = 127;

j=[0:1:N];

x = cos((j+0.5)*pi/(N+1));

k = [0:1:2*N];
xj = cos((k+0.5)*pi/(N+1));

%Function definition
%fn = 1./(1 + 16*x.^2);%Even Function

fn = x.^3 - x;%Odd Function

%Performing the discrete chebyshev transform
c = ones(N+1,1);
c(1,1) = 2;
c(N+1,1) = 2;

Cmat = zeros(N+1,N+1);

for i=1:length(Cmat(:,1))
    for j=1:length(Cmat(:,1))
        Cmat(i,j) = cos((i-1)*(j-1)*pi/N)/c(j,1);
    end
end

u_k = Cmat*fn';
ck = 2./(c*N);
u_k = ck.*u_k;


lincomb = 0.0;
lincombL = 0.0;
lincombLm = 0.0;

%Checking whether the function is odd or even. First kind chebyshev
%polynomial is used when 'm' is even and Second Kind is used when 'm' is
%odd - this gives us four cases 

if mod(m,2) == 0 && fn(1,N) == fn(1,2)          % m = even, function = even
    chebCoeffs = u_k(1:2:length(u_k));
    lincomb = ChebyshevT_M_even_fn_even(chebCoeffs,xj);
    leg2Coeffs = chebyshevT_Legendre2_M_even_An2plus(length(chebCoeffs)-1)*chebCoeffs;
    lincombL = legendre2Interpolation_M_even_fn_even(leg2Coeffs,xj);
    legMCoeffs = (legendre2_LegendreM_M_even_fn_even(length(leg2Coeffs),m))'*leg2Coeffs;
    lincombLm = legendreMInterpolation_M_even_fn_even(legMCoeffs,m,xj);
    %plot(x,fn,'o',xj,lincomb,'--',xj,lincombL,'-s',xj,lincombLm,':p');
    
elseif mod(m,2) == 0 && fn(1,N) == -fn(1,2)     %m = even, function = odd
    chebCoeffs = u_k(2:2:N);
    lincomb = ChebyshevT_M_even_fn_odd(chebCoeffs,xj);
    leg2Coeffs = chebyshevT_Legendre2_M_even_An2minus(length(chebCoeffs))*chebCoeffs;
    lincombL = legendre2Interpolation_M_even_fn_odd(leg2Coeffs,xj);
    legMCoeffs = (legendre2_LegendreM_M_even_fn_odd(length(leg2Coeffs),m))'*leg2Coeffs;
    lincombLm = legendreMInterpolation_M_even_fn_odd(legMCoeffs,m,xj);
    %plot(x,fn,'o',x,lincomb,'--',x,lincombL,':s',x,lincombLm,':p');
    
elseif mod(m,2) ~= 0 && fn(1,N) == fn(1,2)      %m = odd, function = even
    chebCoeffs = u_k(1:2:length(u_k));
    lincomb = ChebyshevU_M_odd_fn_even(chebCoeffs,xj);
    leg1Coeffs = chebyshevU_Legendre1_M_odd_An1plus(length(chebCoeffs))*chebCoeffs;
    lincombL = legendre1Interpolation_M_odd_fn_even(leg1Coeffs,xj);
    legMCoeffs = legendre1_LegendreM_M_odd_fn_even(length(leg1Coeffs),m)'*leg1Coeffs;
    lincombLm = legendreMInterpolation_M_odd_fn_even(legMCoeffs,m,xj);
    %plot(x,fn,'o',x,lincomb,'--',x,lincombL,':s',x,lincombLm,':p');
    
else                                            %m = odd, function = odd
    chebCoeffs = u_k(2:2:N);
    lincomb = ChebyshevU_M_odd_fn_odd(chebCoeffs,xj);
    leg1Coeffs = chebyshevU_Legendre1_M_odd_An1minus(length(chebCoeffs))*chebCoeffs;
    lincombL = legendre1Interpolation_M_odd_fn_odd(leg1Coeffs,xj);
    legMCoeffs = legendre1_LegendreM_M_odd_fn_odd(length(leg1Coeffs),m)'*leg1Coeffs;
    lincombLm = legendreMInterpolation_M_odd_fn_odd(legMCoeffs,m,xj);
    %plot(x,fn,'o',x,lincomb,'--',x,lincombL,':s',x,lincombLm,':p');
end
figure(1)
plot(x,fn,'o');
hold on;
plot(xj,lincomb,'--','LineWidth',2);
hold on;
plot(xj,lincombL,'g-s','MarkerSize',5);
hold on;
plot(xj,lincombLm,':','LineWidth',2);
legend('Original Function','Chebyshev Interpolation',...
    'Lower order Legendre Interpolation','m_{th} order Legendre Interpolation');