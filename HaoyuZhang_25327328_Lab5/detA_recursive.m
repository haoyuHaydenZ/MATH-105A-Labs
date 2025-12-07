function d = detA_recursive(A)
%   d = detA_recursive(A) returns the determinant of a square matrix A
%   using the recursive definition (Laplace expansion).
%   Base cases: n = 1 and n = 2 handled explicitly.
%   For n >= 3, expands along a row/column with the most zeros to reduce work.

    if ndims(A) ~= 2 || size(A,1) ~= size(A,2)
        error('detA_recursive:InputNotSquare', 'Input must be a square matrix.');
    end
    n = size(A,1);

    % ---- base cases ----
    if n == 0
        d = 1;  % conventionally det([]) = 1
        return
    elseif n == 1
        d = A(1,1);
        return
    elseif n == 2
        d = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        return
    end

    if any(all(A == 0, 2)) || any(all(A == 0, 1))
        d = 0;
        return
    end

    zeroRows = sum(A == 0, 2);
    zeroCols = sum(A == 0, 1);
    [zr, rIdx] = max(zeroRows);       
    [zc, cIdx] = max(zeroCols);       

    d = 0;

    if zr >= zc
        i = rIdx;
        for j = 1:n
            aij = A(i,j);
            if aij ~= 0 % (-1)^(i+j) * a_ij * det(M_ij)
                cofactorSign = (-1)^(i+j);
                Mij = A([1:i-1,i+1:n], [1:j-1,j+1:n]);
                d = d + cofactorSign * aij * detA_recursive(Mij);
            end
        end
    else
        j = cIdx;
        for i = 1:n
            aij = A(i,j);
            if aij ~= 0
                cofactorSign = (-1)^(i+j);
                Mij = A([1:i-1,i+1:n], [1:j-1,j+1:n]);
                d = d + cofactorSign * aij * detA_recursive(Mij);
            end
        end
    end
end
