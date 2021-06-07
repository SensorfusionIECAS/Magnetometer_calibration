% normalize 3*3 rotation matrix
function A = matrix_normal(B)
A = zeros(3,3);
A(:,1) = B(:,1) / norm(B(:,1));
A(:,2) = B(:,2) / norm(B(:,2));
A(:,3) = B(:,3) / norm(B(:,3));