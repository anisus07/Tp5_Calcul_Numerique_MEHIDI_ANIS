function [x]=Richardson(A,b,x0,k)
n=size(A,"r")
x=x0
for i=1:k
    for j=1:n
        x(j) = ((b(j)-A(j,1:n)* x )* A(j,j)) + x(j);
    end
end
endfunction

function [x, erreur] = test_Richardson()
    
    A = [2 -1 0; -1 2 -1;0 -1 2]
    b = [5;6;2]
    x0 = zeros(3,1)
    k = 50
    [x] = Richardson(A,b,x0,k)
    erreur = norm((A*x-b)/x)

endfunction
