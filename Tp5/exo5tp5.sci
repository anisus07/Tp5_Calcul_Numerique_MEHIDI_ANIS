function [L,U] =factorisation(A)
    [m,n]=size(A);
    
    for i= 1:n
        L(i,1) = A(i,1);
    end
    
    for j=1:m
        U(1,j)=A(1,j)/L(1,1);
    end
    
    for j=2:n
        for i=j:n
            r=0.0;
            
            for k= 1:(j-1)
                r=r+L(i,k)*U(k,j);
                
            end
            L(i,j)=A(i,j)-r 
        end
     U(j,j)=1;
     for i=(j+1):n
         r=0.0;
       for k=1:(j-1)
           r=r+L(j,k)*U(k,i);
       end
       U(j,i)=(A(j,i)-r)/L(j,j);                
    end
    end
endfunction
t = zeros(10);
x=zeros(10);
d = [10:10:100];
    for n = d
    
        i = n / 10;
    
        // Init matrix
        A = rand(n, n);
        
        disp(A);
        x(i)=n*n;
        tic;
        [L,U] =factorisation(A);
        t(i) = toc();
        
    end

   xtitle("Factorisation LU", "n", "time");
plot(d, [t x]);
legend(["Factorisation LU" "Complexité"], 2);
xs2png(0, "./graphe/facto.png");
clf();
