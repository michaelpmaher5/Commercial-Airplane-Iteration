function solutions = solveForWto(A,B,C,D)
    %%syms x
    %%solutions = vpasolve(A*x^1.195 + B*x^0.235 + C*x + D == 0, x);
    %%k = find(values > 100000 & values < 400000,1);
    %%solutions = values(k);
    x = 50000:50:500000;
    %%checkmatrix = zeros(size(x));
    test = abs(A*x.^1.195 + B*x.^0.235 + C*x + D);
    %%solutions = 500;
    [lookfor,index] = min(test,[],'omitnan');
    solutions = x(index);
    %%for i = 50000:500:900000
       %% test2 = abs(A*i^1.195 + B*i^0.235 + C*i + D);
        %%if(abs(test2 - lookfor) < 10)
          %%  solutions = i;
       %% end
   %% end
    
    
   %% for x = 50000:500:900000
     %%   test = A*x^1.195 + B*x^0.235 + C*x + D;
       %% checkmatrix(x,2) = abs(test);
        %%if(abs(test) < 1000)
           %% solutions = x;
       %% end
        
    %%end
    
    
end