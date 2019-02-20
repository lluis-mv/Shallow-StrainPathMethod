
function [] = PostProcessResults(XFILE, matFILE)

global IAmAMechanicalEngineer
IAmAMechanicalEngineer = true;

if (nargin ==  0)
    XFILE = 'StrainPath';
    matFILE = 'StrainPath.mat';
end

global D_over_T

%Not the perfect way of doing it.
load(matFILE)


minStrain = -0.05;
maxStrain = - minStrain;

rr = rr(:,1:jj); zz = -zz(:,1:jj);
disp_radial = disp_radial(:,1:jj); disp_vertical = -disp_vertical(:,1:jj);
ZZ = -ZZ(:,1:jj);
RR = RR(:,1:jj);
eTheta = eTheta(:, 1:jj);
eZ = eZ(:, 1:jj);
eRZ = eRZ(:, 1:jj);
eR = eR(:, 1:jj);
vr = vr(:,1:jj); vz = -vz(:,1:jj);
FlagShallow = exist('vr1', 'var');

if ( FlagShallow)
    vr1 = vr1(:,1:jj); vz1 = -vz1(:,1:jj);
    vr2 = vr2(:,1:jj); vz2 = -vz2(:,1:jj);
    vr3 = vr3(:,1:jj); vz3 = -vz3(:,1:jj);
end

U = ConvertToVector(disp_radial, disp_vertical);
[X, T] = ConvertToFEMMesh( rr, zz, U);
V = ConvertToVector(vr, vz);

if ( FlagShallow )
    V1 = ConvertToVector(vr1, vz1);
    V2 = ConvertToVector(vr2, vz2);
    V3 = ConvertToVector(vr3, vz3);
end

if ( FlagShallow == false )
    % Scale The Model
    Rpfem = 0.0375;
    X = X*Rpfem;
    U = U*Rpfem;
    V = V*Rpfem;
    X(:,2) = X(:,2)-3*Rpfem;
else
    % Scale The Model
    Rpfem = 0.0375;
    X = X*Rpfem;
    U = U*Rpfem;
    V = V*Rpfem;
    V1 = V1*Rpfem;
    V2 = V2*Rpfem;
    V3 = V3*Rpfem;
end

if ( FlagShallow )
    WriteToGid(X, T, U, V, XFILE, V1, V2, V3);
else
    WriteToGid(X, T, U, V, XFILE);
end


function [] = WriteToGid( X, T, U, V, XFILE, V1, V2, V3)

nNodes = size(X,1);
nElem = size(T, 1);
% write The Mesh
fid = fopen( ['a', XFILE, '.msh'], 'w' );

fprintf(fid, 'MESH  dimension 2  ElemType Quadrilateral Nnode 4 \n');
fprintf(fid, 'Coordinates \n ');
for i = 1:nNodes
    fprintf(fid, '%i %e %e %e \n', [i, X(i,:), 0]);
end
fprintf(fid, 'END Coordinates \n');


fprintf(fid, 'elements \n');
for i = 1:nElem
    fprintf( fid, '%i %i %i %i %i \n', [i, T(i,:)]);
end
fprintf(fid, 'END elements \n');
fclose(fid);

fid = fopen(['a', XFILE, '.res'], 'w' );
fprintf(fid, 'Gid Post Results File 1.0 \n');



fprintf( fid, 'GaussPoints "GP"  ElemType Quadrilateral  \n');
fprintf(fid, ' Number Of Gauss Points: 4 \n');
fprintf(fid, 'Natural Coordinates: Internal \n');
fprintf(fid, 'End GaussPoints \n');

time = num2str(rand());

fprintf( fid, ['Result "DISPLACEMENT" "Kratos" ', time,' Vector OnNodes \n']);
fprintf( fid, 'ComponentNames "X-DISPLACEMENT", "Y-DISPLACEMENT", "Z-DISPLACEMENT" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, U(i,:), 0] );
end
fprintf( fid, 'End Values \n');


fprintf( fid, ['Result "Velocity" "Kratos" ', time,'  Vector OnNodes \n']);
fprintf( fid, 'ComponentNames "X-Velocity", "Y-Velocity", "Z-Velocity" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, V(i,:), 0] );
end
fprintf( fid, 'End Values \n');

if ( nargin == 8)
        fprintf( fid, ['Result "Velocity1" "Kratos" ', time,'  Vector OnNodes \n']);
    fprintf( fid, 'ComponentNames "X-Velocity1", "Y-Velocity1", "Z-Velocity1" \n');
    fprintf( fid, 'Values \n');
    for i = 1:nNodes
        fprintf( fid, '%i %e %e %e \n', [i, V1(i,:), 0] );
    end
    fprintf( fid, 'End Values \n');
    fprintf( fid, ['Result "Velocity2" "Kratos" ', time,'  Vector OnNodes \n']);
    fprintf( fid, 'ComponentNames "X-Velocity2", "Y-Velocity2", "Z-Velocity2" \n');
    fprintf( fid, 'Values \n');
    for i = 1:nNodes
        fprintf( fid, '%i %e %e %e \n', [i, V2(i,:), 0] );
    end
    fprintf( fid, 'End Values \n');
    fprintf( fid, ['Result "Velocity3" "Kratos" ', time,'  Vector OnNodes \n']);
    fprintf( fid, 'ComponentNames "X-Velocity3", "Y-Velocity3", "Z-Velocity3" \n');
    fprintf( fid, 'Values \n');
    for i = 1:nNodes
        fprintf( fid, '%i %e %e %e \n', [i, V3(i,:), 0] );
    end
    fprintf( fid, 'End Values \n');
end

fprintf(fid, ['Result "Small_strain_tensor" "Kratos" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeStrains( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n');



fprintf(fid, ['Result "ALMANSI_STRAIN_TENSOR" "Kratos" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeAlmansi( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n ');

fprintf(fid, ['Result "difference_strain_tensor" "Kratos" ', time,'  Matrix OnGaussPoints "GP" \n ']);
fprintf(fid, ' values \n ');
for i = 1:nElem
    Celem = T(i,:);
    [S1, S2, S3, S4] = ComputeDifference( X(Celem, :), U(Celem,:));
    fprintf( fid, '%i %e %e %e %e %e %e \n', [i, S1] );
    fprintf( fid, '   %e %e %e %e %e %e \n', S2 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S3 );
    fprintf( fid, '   %e %e %e %e %e %e \n', S4 );
end
fprintf(fid, ' END values \n ');


fclose(fid);

system(['mv a', XFILE, '.res ',  XFILE, '.res ']);
system(['mv a', XFILE, '.msh ',  XFILE, '.msh ']);

function [S1, S2, S3, S4] = ComputeDifference( xe, ue)


A = StrainsDifference( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


function [S1, S2, S3, S4] = ComputeAlmansi( xe, ue)

global IAmAMechanicalEngineer
sign = 1;

if (IAmAMechanicalEngineer)
    sign = - sign;
end

A = StrainsAlmansi( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = sign*A(:,1)';
S2 = sign*A(:,2)';
S3 = sign*A(:,3)';
S4 = sign*A(:,4)';


function [S1, S2, S3, S4] = ComputeStrains( xe, ue)
       
global IAmAMechanicalEngineer
sign = 1;

if (IAmAMechanicalEngineer)
    sign = - sign;
end

A = StrainsEpsilon( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = sign*A(:,1)';
S2 = sign*A(:,2)';
S3 = sign*A(:,3)';
S4 = sign*A(:,4)';


function []=PlotTube()
global D_over_T
D = 2;
t = D/D_over_T;
hold on
x = [1, 1-t, 1-t, 1];
y = [0, 0 , 4, 4];
z = 100*[1,1,1,1];
patch(x,y,z, 'k')


x=[]; y=[]; z=[];

for i = 1:100
    a = i/100*2*pi;
    x = [x, (1-t/2)+t/2*sin(a)];
    y = [y, t/2*cos(a)];
    z = [100, z];
end
patch(x,y,z,'k')
hold off



function [X, T ] = ConvertToFEMMesh(rr, zz, U)

global D_over_T

nNodes = size(rr, 1)*size(rr, 2);
X = zeros(nNodes, 2);


for i = 1:size(rr,1)
    for j = 1:size(rr,2)
        iNod = i + (j-1)*size(rr,1);
        X(iNod, 1) = rr(i,j);
        X(iNod, 2) = zz(i,j);
    end
end

T = [];

for i = 1:size(rr,1)-1
    for j = 1:size(rr,2)-1
        nNod1 = i + (j-1)*size(rr,1);
        nNod2 = i +(j)*size(rr,1);
        nNod3 = i +1 +(j)*size(rr,1);
        nNod4 = i +1 +(j-1)*size(rr,1);
        
        nNod = [nNod1, nNod2, nNod3, nNod4, nNod1]';
        
        xNod = X(nNod, :);
        
        XNod =  X(nNod, :) - U(nNod,:);
        
        a = polyarea(xNod(:,1), xNod(:,2));
        A = polyarea(XNod(:,1), XNod(:,2));
        
        ratio = abs( (A-a) / A);
        
        
        if (  ratio > 0.25 )
            % do nothing
        else
            T = [T;
                nNod1, nNod2, nNod3, nNod4];
        end
    end
end


function U  = ConvertToVector( u1, u2)

nNodes = size(u1, 1)*size(u1, 2);
U = zeros(nNodes, 2);

for i = 1:size(u1,1)
    for j = 1:size(u1,2)
        nNod = i + (j-1)*size(u1,1);
        U(nNod, 1) = u1(i,j);
        U(nNod, 2) = u2(i,j);
    end
end