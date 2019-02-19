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
D_over_T = 41;


disp('++Generate ODE source term::DONE')
disp('Generate strain files')
tic
% SymbolicComputeStrains();
toc
disp('++Generate strain files::DONE')



[RR, ZZ] = meshgrid([linspace(1e-6, 1, 20), linspace(1.05, 2.5, 60)],  ...
    linspace(0,6.5,50) );
% allocatate memory
rr = 0*RR; zz = rr; disp_radial = rr; disp_vertical = rr;
eR = rr; eZ = rr; eTheta = rr; eRZ = rr;
vr = rr; vz = rr;

InitialZOfTube = 0;
FinalZOfTube = 5;

waitBar = waitbar(0, 'Computing displacements and strains');



for jj = 1:size(RR,2)
    for ii = 1:size(RR,1)
        R = RR(ii,jj); Z = ZZ(ii,jj);
        
        [r, z] =IntegrateDisplacements( InitialZOfTube, ...
            FinalZOfTube, R, Z);
        
        rr(ii,jj) = r; 
        zz(ii,jj) = z;
        [vr(ii,jj), vz(ii,jj)] = EvaluateVelocity(r, z, FinalZOfTube);
    end
    
    waitbar(jj/size(RR,2),waitBar, 'Computing displacements and strains');
    pause(0.0001)
    disp_radial = rr - RR;
    disp_vertical = zz - ZZ;
    
    save('variablesT.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
        'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
    if ( jj > 3)
        PostProcessResults();
        disp('PostProcess')
    end
end


save('variables.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
    'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
close(waitBar)


function [r, z] = IntegrateDisplacements(hIni, hEnd, R, Z)

initialCondition = zeros(2,1);
initialCondition(1) = R;
initialCondition(2) = Z;

SourceFunction = @(t,x) AuxiliarFunction(t,x);
options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
[t, xx] = ode45( SourceFunction , [hIni, hEnd], initialCondition, options);


r = xx(end,1);
z = xx(end,2);
    


function [vr, vz] = EvaluateVelocity(r, z, FinalZofTube)

[a] = SourceTermShallowStrainPath(FinalZofTube, r, z);
vr = a(1);
vz = a(2);

function [dxdt] = AuxiliarFunction(t, x)
r = x(1);
z = x(2);
[v] = SourceTermShallowStrainPath(t,r,z);

dxdt = v;
