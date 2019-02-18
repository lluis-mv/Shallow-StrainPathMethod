%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab symbolic file to generate the stiffness matrix 
% of an axisymmetric potential flow (in terms of the potential)
%      LMonforte - February 2019
%

% lets assume 2D
ndim = 2;
nNodes = 3;

% parametric coordinates 
syms alfa real; 
syms beta real;

%
U = sym( 'u', [nNodes 1], 'real');
X = sym( 'x', [3,2], 'real');



%Shape funcitons
Nsmall = sym( 'n', [nNodes 1], 'real');
Nsmall(1) = 1 - alfa - beta;
Nsmall(2) = alfa; 
Nsmall(3) = beta;

%Derivative shape functions with respect to parametric
Nsmall_chi = [-1 -1; 1 0; 0 1];


% Continuous solution and its gradient in terms of X
Ucontinuous = Nsmall'*U;

J = sym(zeros(2,2));
for j = 1:nNodes
    J = J + [ Nsmall_chi(j,1) * X(j,1), Nsmall_chi(j,1) * X(j,2);
        Nsmall_chi(j,2) * X(j,1), Nsmall_chi(j,2) * X(j,2)];
end

dN_dX = J \ Nsmall_chi';
dU_dX = dN_dX * U;

Guess = sym('gg', [nNodes 1], 'real');
dGuess_dX = dN_dX*Guess;


otherU = [diff(Ucontinuous, alfa); diff(Ucontinuous, beta) ];
otherU = J\otherU ;


% Define the weak form and its tangent matrix
r = Nsmall'*X(:,1);
weakForm = (dGuess_dX' * dU_dX) * det(J)*2*pi*r;
weakForm = int(weakForm, beta, 0, 1-alfa);
weakForm = int(weakForm, alfa, 0, 1);

residual = sym(zeros(nNodes, 1));
for i = 1:nNodes
    residual(i) = diff(weakForm, Guess(i));
end

K = sym(zeros(nNodes,nNodes));
for i = 1:nNodes
    for j = 1:nNodes
        K(i,j) = diff(residual(i), U(j));
    end
end

velocity = dU_dX';


return;

f = matlabFunction(residual, K, velocity, 'File','StrainPathProblem', 'Vars', [X, U]);

% Second part. Try To Compute the deformation Gradient F;
XRef = sym( 'xRef', [3,2], 'real');
J = sym(zeros(2,2));
for j = 1:nNodes
    J = J + [ Nsmall_chi(j,1) * XRef(j,1), Nsmall_chi(j,1) * XRef(j,2);
        Nsmall_chi(j,2) * XRef(j,1), Nsmall_chi(j,2) * XRef(j,2)];
end

dN_dXRef = inv(J) * Nsmall_chi';
r = sum(X(:,1))/3;
R = sum(XRef(:,1))/3;

F = sym(zeros(3,3));
for i = 1:2
    for j = 1:2
        for nn = 1:3
            F(i,j) = F(i,j) + X(nn,i) * dN_dXRef(j,nn);
        end
    end
end
F(3,3) = r/R;
F = simplify(F);
detF = det(F);
detF = simplify(detF);


f = matlabFunction(F,detF, 'File','StrainPathComputeF', 'Vars', [XRef, X]);