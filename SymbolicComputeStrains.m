%=========================================================================
%  Symbolic Matlab to compute the deformation and strains once the
%  displacement are known
%
%  Following FEM methodology
%
%
%  Barcelona, 13 February 2019
%=========================================================================
function [] = SymbolicComputeStrains()

% They are defined here but can not be changed
nDim = 2;
nNodes = 4;

syms chi real;
syms eta real;

U = sym('u', [nNodes, nDim], 'real');
X = sym('x', [nNodes, nDim], 'real');

% Shape functions
N = sym(zeros( nNodes, 1));
N(1) = (1-chi)*(1-eta);
N(2) = (1+chi)*(1-eta);
N(3) = (1+chi)*(1+eta);
N(4) = (1-chi)*(1+eta);
N = N/4;

Ucontinuous = N'*U;
r = N' * X(:,1);

dNdChi = [diff(N, chi), diff(N, eta)];

dXdChi = sym(zeros(nDim, nDim));
for j = 1:nDim
    for k = 1:nDim
        for i = 1:nNodes
            dXdChi(j,k) = dXdChi(j,k) + X(i,j)*dNdChi(i,k);
        end
    end
end

dNdX = dNdChi* inv(dXdChi);

dUdX = sym(zeros(nDim, nDim));
for j = 1:nDim
    for k = 1:nDim
        for i = 1:nNodes
            dUdX(j,k) = dUdX(j,k) + U(i,j)*dNdX(i,k);
        end
    end
end

epsilon = 1/2*( dUdX  + dUdX');
epsilon = [epsilon, zeros(2,1);
    zeros(1,2), Ucontinuous(1)/ r];

eVoigt = CompressInVoigtNotation( epsilon);

[eP1, eP2, eP3, eP4] = EvaluateStrainAtGP( eVoigt);

matlabFunction( [eP1, eP2, eP3, eP4], 'File', 'StrainsEpsilon', 'Vars', [X, U]);

AlmansiTensor = 1/2*( dUdX + dUdX');

AlmansiTensor(1,1) = AlmansiTensor(1,1) - 0.5*( dUdX(1,1)^2 + dUdX(1,2)^2);
AlmansiTensor(2,2) = AlmansiTensor(2,2) - 0.5*( dUdX(2,1)^2 + dUdX(2,2)^2);

AlmansiTensor(1,2) = AlmansiTensor(1,2) - 0.5* ( dUdX(2,1)*dUdX(1,1) + dUdX(1,2)*dUdX(2,2));
AlmansiTensor(2,1) = AlmansiTensor(1,2);

AlmansiTensor = [AlmansiTensor, zeros(2,1);
    zeros(1,3)];
AlmansiTensor(3,3) = Ucontinuous(1)/r - 1/2* ( Ucontinuous(1)/r)^2;

aVoigt = CompressInVoigtNotation( AlmansiTensor);

[aP1, aP2, aP3, aP4] = EvaluateStrainAtGP( aVoigt);

matlabFunction( [aP1, aP2, aP3, aP4], 'File', 'StrainsAlmansi', 'Vars', [X, U]);

matlabFunction( [eP1-aP1, eP2-aP2, eP3-aP3, eP4-aP4], 'File', 'StrainsDifference', 'Vars', [X, U]);

function eVoigt = CompressInVoigtNotation( epsilon)



eVoigt = sym(zeros(6,1));

for i = 1:3
    eVoigt(i) = -epsilon(i,i);
end
eVoigt(4) = epsilon(1,2)+epsilon(2,1);
eVoigt(5) = epsilon(1,3)+epsilon(3,1);
eVoigt(6) = epsilon(2,3)+epsilon(3,2);


function [eP1, eP2, eP3, eP4] = EvaluateStrainAtGP( Strain)
syms chi real;
syms eta real;

number = sqrt(3)/3;

eP1 = subs( Strain, chi, -number);
eP1 = subs( eP1, eta, -number);

eP2 = subs( Strain, chi, number);
eP2 = subs( eP2, eta, -number);

eP3 = subs( Strain, chi, number);
eP3 = subs( eP3, eta, number);

eP4 = subs( Strain, chi, -number);
eP4 = subs( eP4, eta, number);