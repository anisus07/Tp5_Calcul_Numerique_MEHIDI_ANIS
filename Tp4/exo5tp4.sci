function [y]=matrice_vecteur_csr(col,val,row,x)

    n = size(row,2)-1
    y = zeros(n,1)

for i = 1 : n
  for j = row(i) : row(i+1)-1
    y(i) = y(i)+ val(j)*x(col(j));
  end
end

endfunction
