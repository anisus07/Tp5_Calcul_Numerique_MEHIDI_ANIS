function [x,e]=jacobi(A,b,n,epsilon)

     N= -(triu(A,1)+tril(A,-1))
     M= diag(diag(A))

    x=zeros(n,1);
    e=norm(A*x-b)


    while(e>epsilon)
     x = inv(M)*(N*x + b)
     e=norm(A*x-b)
    end

endfunction


function [x,e]=gauss_seidel(A,b,n,epsilon)
    
    z=0
    x=zeros(n,1)
    
    e=1
    while (e>epsilon)
        for i=1:n
            y=0
            for  j=1:i-1
               y=y+A(i,j)*x(j)
            end
            z=0
            for j=i+1:n
                z=z+A(i,j)*x(j)
            end 
            x(i)=(1/A(i,i))*(b(i)-z-y)
        end
        e=norm(b-A*x)
        
    end
        
    
    
endfunction

