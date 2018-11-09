clc
clear 

J = 4; 

% For Dirichlet:

% e = ones((J-1)*(J-1), 1);
% A = spdiags([-e -e 4*e -e -e], [-(J-1), -1, 0, 1, J-1], (J-1)*(J-1), (J-1)*(J-1));
% B = full(A)
% 
% for i = 1:(J-1)-1
%    A( i*(J-1)+1, i*(J-1) ) = 0;
%    A( i*(J-1), i*(J-1)+1 ) = 0;
% end
% B2 = full(A)


% For periodic:
e = ones(J*J, 1);
A = spdiags([-e -e -e 4*e -e -e -e], [-(J-1)*J ,-J, -1, 0, 1, J, (J-1)*J], J*J, J*J);
B = full(A)

for i = 1:J-1
   A( (i*J)+1, i*J ) = 0;
   A( i*J, (i*J)+1 ) = 0;
   
   A( i*J, (i-1)*J + 1 ) = -1;
   A( (i-1)*J + 1, i*J ) = -1;
   
end
A( J*J, (J-1)*J + 1 ) = -1;
A( (J-1)*J + 1, J*J ) = -1;

B2 = full(A)