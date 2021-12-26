function [x,e]=jacobi(A,b,n,epsilon)

    M=-(triu(A,1)+tril(A,-1))
    temp=diag(A)
    N=zeros(n,n)
    for i=1:n
        N(i,i)=temp(i)
    end
    
    x=zeros(n,1);
    B=inv(N)*M
    v=inv(N)*b
    temp=zeros(n,1)
    ev=zeros(n,1)
    e=1
   
    
    while(e>epsilon)
      ev=zeros(n,1)
      temp=zeros(n,1)
      for i=1:n
           for j=1:n
           temp(i)=temp(i)+B(i,j)*x(j)
           end 
      end
      temp(i)=temp(i)+v(i)
    e=norm(b-A*x)   
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

