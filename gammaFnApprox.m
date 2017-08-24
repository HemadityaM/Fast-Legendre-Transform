function [val1,val2] = gammaFnApprox(temp1,temp2)
        %temp1 = k-l-1;
        %temp2 = (2*k+2*l+1)/2;
        G1 = (gamma(temp1+0.5))/(gamma(temp1+1));
        G2 = (gamma(0.5+temp2))/(gamma(1+temp2));
            if isnan(G1)
                %fprintf('l+1 = %d, k+1 = %d \n',l+1,k+1);
                val1 = order6SeriesExpansion(temp1);
            else
                val1 = G1;
            end
            if isnan(G2)
                %fprintf('l+1 = %d, k+1 = %d \n',l+1,k+1);
                val2 = order6SeriesExpansion(temp2);
            else
                val2 = G2;
            end
end