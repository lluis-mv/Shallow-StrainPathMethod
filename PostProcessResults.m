
function [] = PostProcessResults()
clear all;

global D_over_T

%Not the perfect way of doing it.
load('variablesT.mat')


minStrain = -0.05;
maxStrain = - minStrain;

rr = rr(:,1:jj);
zz = -zz(:,1:jj);
disp_radial = disp_radial(:,1:jj);
disp_vertical = -disp_vertical(:,1:jj);
ZZ = -ZZ(:,1:jj);
RR = RR(:,1:jj);
eTheta = eTheta(:, 1:jj);
eZ = eZ(:, 1:jj);
eRZ = eRZ(:, 1:jj);
eR = eR(:, 1:jj);
vr = vr(:,1:jj);
vz = -vz(:,1:jj);


[X, T, U, V] = ConvertToFEMMesh( rr, zz, disp_radial, disp_vertical, vr, vz);

WriteToGid(X, T, U, V, 'Results');


function [] = WriteToGid( X, T, U, V, XFILE)

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


fprintf( fid, 'Result "Displacements" "StrainPathMethod" 1 Vector OnNodes \n');
fprintf( fid, 'ComponentNames "X-Displ", "Y-Displ", "Z-Displ" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, U(i,:), 0] );
end
fprintf( fid, 'End Values \n');


fprintf( fid, 'Result "Velocity" "StrainPathMethod" 1 Vector OnNodes \n');
fprintf( fid, 'ComponentNames "X-Velocity", "Y-Velocity", "Z-Velocity" \n');
fprintf( fid, 'Values \n');
for i = 1:nNodes
    fprintf( fid, '%i %e %e %e \n', [i, V(i,:), 0] );
end
fprintf( fid, 'End Values \n');

fprintf(fid, 'Result "Small_strain_tensor" "StrainPathMethod" 1 Matrix OnGaussPoints "GP" \n ');
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



fprintf(fid, 'Result "Almansi_strain_tensor" "StrainPathMethod" 1 Matrix OnGaussPoints "GP" \n ');
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

fprintf(fid, 'Result "difference_strain_tensor" "StrainPathMethod" 1 Matrix OnGaussPoints "GP" \n ');
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


A = StrainsAlmansi( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


function [S1, S2, S3, S4] = ComputeStrains( xe, ue)


A = StrainsEpsilon( xe(1,1), xe(2,1), xe(3,1), xe(4,1), ...
    xe(1,2), xe(2,2), xe(3,2), xe(4,2), ...
    ue(1,1), ue(2,1), ue(3,1), ue(4,1), ...
    ue(1,2), ue(2,2), ue(3,2), ue(4,2));

S1 = A(:,1)';
S2 = A(:,2)';
S3 = A(:,3)';
S4 = A(:,4)';


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



function [X, T, U, V ] = ConvertToFEMMesh(rr, zz, disp_radial, ...
    disp_vertical, vr, vz)

global D_over_T

nNodes = size(rr, 1)*size(rr, 2);
X = zeros(nNodes, 2);
U = zeros(nNodes, 2);
V = zeros(nNodes, 2);

for i = 1:size(rr,1)
    for j = 1:size(rr,2)
        nNod = i + (j-1)*size(rr,1);
        X(nNod, 1) = rr(i,j);
        X(nNod, 2) = zz(i,j);
        U(nNod, 1) = disp_radial(i,j);
        U(nNod, 2) = disp_vertical(i,j);
        V(nNod, 1) = vr(i,j);
        V(nNod, 2) = vz(i,j);
    end
end

T = [];

for i = 1:size(rr,1)-1
    for j = 1:size(rr,2)-1
        nNod1 = i + (j-1)*size(rr,1);
        nNod2 = i +(j)*size(rr,1);
        nNod3 = i +1 +(j)*size(rr,1);
        nNod4 = i +1 +(j-1)*size(rr,1);
        xNod1 = X(nNod1, :);
        xNod2 = X(nNod2, :);
        xNod3 = X(nNod3, :);
        xNod4 = X(nNod4, :);
        meanN = 1/4*(xNod1 + xNod2 + xNod3 + xNod4);
%         if (  (meanN(2) > 0) && ( meanN(1) < 0.99 && meanN(1) > 1-0.98*2/D_over_T) )
%             % do nothing
%         else
            T = [T;
                nNod1, nNod2, nNod3, nNod4];
%         end
    end
end


