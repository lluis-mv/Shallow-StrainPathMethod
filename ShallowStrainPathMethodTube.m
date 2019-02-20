%=========================================================================
%  Matlab to compute the velocities and strains Strains due to   
%  the penetration of a tube using the Strain Path Method (Baligh, 1985)  
%
%  Following the approach of Sagaseta et al (1991)
%
%  This file is responsible of the numerical integration whereas the other
%  file defines the velocity and temporal derivative of strains
%  simbolically and write them into a file 
%  
%  Shallow strain path method!
%
%  Barcelona, 13 February 2019
%=========================================================================

% V0: To postprocess, run the "SecondPart". It might be run during the
% computation of the first part (in another matlab)

function []=ShallowStrainPathMethodTube()

% Define the geometry of the tube
%   the external radius, velocity of the tube are set to one
global D_over_T
D_over_T = 10;


disp('++Generate ODE source term::DONE')
disp('Generate strain files')
tic
% SymbolicComputeStrains();
toc
disp('++Generate strain files::DONE')



[RR, ZZ] = meshgrid([linspace(1e-6, 1, 15), linspace(1.05, 4, 20)],  ...
    linspace(0.1,6,35) );
% allocatate memory
rr = 0*RR; zz = rr; disp_radial = rr; disp_vertical = rr;
eR = rr; eZ = rr; eTheta = rr; eRZ = rr;
vr = rr; vz = rr; vr1 = rr; vz1 = rr;
vr2 = rr; vz2 = rr; vr3 = rr; vz3 = rr;

InitialZOfTube = 0;
FinalZOfTube = 3;

% waitBar = waitbar(0, 'Computing displacements and strains');



for jj = 1:size(RR,2)
    tic
    for ii = 1:size(RR,1)
        R = RR(ii,jj); Z = ZZ(ii,jj);
        
        [r, z] =IntegrateDisplacements( InitialZOfTube, ...
            FinalZOfTube, R, Z);
        
        rr(ii,jj) = r; 
        zz(ii,jj) = z;
        [vr(ii,jj), vz(ii,jj)] = EvaluateVelocity(r, z, FinalZOfTube);
        [vr(ii,jj), vz(ii,jj), vr1(ii,jj), vz1(ii,jj), vr2(ii,jj), vz2(ii,jj), ...
            vr3(ii,jj), vz3(ii,jj) ] = EvaluateVelocity(r, z, FinalZOfTube);
    end
    
%     waitbar(jj/size(RR,2),waitBar, 'Computing displacements and strains');
    pause(0.0001)
    disp_radial = rr - RR;
    disp_vertical = zz - ZZ;
    
    save('Shallow.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
        'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
        'vr1', 'vz1', 'vr2', 'vz2', 'vr3', 'vz3', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
    if ( jj > 3)
        PostProcessResults('ShallowStrainPath', 'Shallow.mat');
        disp(['PostProcess ', num2str(jj), ' of ', num2str( size(RR,2) )] )
        toc;
    end
end


save('Shallow.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
    'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
    'vr1', 'vz1', 'vr2', 'vz2', 'vr3', 'vz3', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
% close(waitBar)


function [r, z] = IntegrateDisplacements(hIni, hEnd, R, Z)

initialCondition = zeros(2,1);
initialCondition(1) = R;
initialCondition(2) = Z;

SourceFunction = @(t,x) AuxiliarFunction(t,x);
options = odeset('RelTol',1e-6, 'AbsTol', 1e-6);
[t, xx] = ode45( SourceFunction , [hIni, hEnd], initialCondition, options);


r = xx(end,1);
z = xx(end,2);
    


function [vr, vz, vr1, vz1, vr2, vz2, vr3, vz3] = EvaluateVelocity(r, z, FinalZofTube)
if (nargout == 2)
    [a] = SourceTermShallowStrainPath(FinalZofTube, r, z);
    vr = a(1);
    vz = a(2);
elseif (nargout == 8)
    [a, b, c, d] = SourceTermShallowStrainPath(FinalZofTube, r, z);
    vr = a(1); vz = a(2);
    vr1 = b(1); vz1 = b(2);
    vr2 = c(1); vz2 = c(2);
    vr3 = d(1); vz3 = d(2);
end

function [dxdt] = AuxiliarFunction(t, x)
r = x(1);
z = x(2);
[v] = SourceTermShallowStrainPath(t,r,z);

dxdt = v;
