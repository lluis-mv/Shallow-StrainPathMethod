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
%
%  Barcelona, 13 February 2019
%=========================================================================

% V0: To postprocess, run the "SecondPart". It might be run during the
% computation of the first part (in another matlab)

function []=StrainPathMethodTube()

% Define the geometry of the tube
%   the external radius, velocity of the tube are set to one
global D_over_T
D_over_T = 10;

% Generate necessary matlab functions
disp('Generate ODE source term')
tic;
SymbolicVelocityDueToTube();
toc
disp('++Generate ODE source term::DONE')
disp('Generate strain files')
tic
SymbolicComputeStrains();
toc
disp('++Generate strain files::DONE')
%********************************************************************
%********************************************************************
%********************** VALIDATION AREA *****************************
%************* compare the CSL with the Closed form *****************
if ( true)
    figure(2105)
    zz = linspace(-20, 20, 1000);
    verticalStrain = -log( 1 + 2 /D_over_T * zz/2 ...
        ./ ( 1 + 4 * ( zz/2).^2).^(3/2));
    plot( verticalStrain , zz,  'k', 'linewidth', 2)
    hold on
    
    [t,XX] =IntegrateDisplacements( -30, 30, 1e-6, 0);
    hold on; plot(XX(:,4), t-XX(:,2), 'linewidth', 1.5)
    ylim([-20,20])
    hold off
    ll = legend('Baligh et al (1987)', 'Numerical integration', ...
        'location', 'best');
    set(ll, 'interpreter', 'latex');
    
    set(gca, 'FontSize', 12)
    grid minor
    drawnow;
    xlabel('Vertical deformation, $\epsilon_z$', 'interpreter', 'latex')
    ylabel('normalized soil position, $z/R$', 'interpreter', 'latex')

end
%****************************************************************
%****************************************************************




[RR, ZZ] = meshgrid([linspace(1e-6, 1, 100), linspace(1.1, 4, 100)],  ...
    linspace(6,-6,100) );
% allocatate memory
rr = 0*RR; zz = rr; disp_radial = rr; disp_vertical = rr;
eR = 0*RR; eZ = rr; eTheta = rr; eRZ = rr;

InitialZOfTube = -25;
FinalZOfTube = 0;

waitBar = waitbar(0, 'Computing displacements and strains');



for jj = 1:size(RR,2)
    for ii = 1:size(RR,1)
        R = RR(ii,jj); Z = ZZ(ii,jj);
        
        [r, z, epsilon] =IntegrateDisplacements( InitialZOfTube, ...
            FinalZOfTube, R, Z);
        
        rr(ii,jj) = r; zz(ii,jj) = z;
        eR(ii,jj) = epsilon(1); eZ(ii,jj) = epsilon(2);
        eTheta(ii,jj) = epsilon(1); eRZ(ii,jj) = epsilon(4);
    end
    
    waitbar(jj/size(RR,2),waitBar, 'Computing displacements and strains');
    pause(0.0001)
    disp_radial = rr - RR;
    disp_vertical = zz - ZZ;
    save('variablesT.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
        'disp_vertical', 'ii', 'jj', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
    if ( jj > 3)
        PostProcessResults();
    end
end


save('variables.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
    'disp_vertical', 'ii', 'jj', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
close(waitBar)


function [r, z, epsilon] = IntegrateDisplacements(hIni, hEnd, R, Z)

initialCondition = zeros(6,1);
initialCondition(1) = R;
initialCondition(2) = Z;

SourceFunction = @(t,x) AuxiliarFunction(t,x);

[t, xx] = ode45( SourceFunction , [hIni, hEnd], initialCondition);

if (nargout == 3)
    r = xx(end,1);
    z = xx(end,2); 
    epsilon = xx(end, 3:6);
elseif (nargout == 2)
    r = t;
    z = xx;
else
    error('not considered case')
end


function [dxdt] = AuxiliarFunction(t, x)
r = x(1);
z = x(2);
[a,b] = SourceTermStrainPath(t,r,z);

dxdt = [a';b'];
