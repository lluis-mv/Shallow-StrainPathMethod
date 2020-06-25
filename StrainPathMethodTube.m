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


function []=StrainPathMethodTube()

% Define the geometry of the tube
%   the external radius, velocity of the tube are set to one

D_over_T = 10;

% Generate necessary matlab functions
disp('Generate ODE source term')
disp('It may take several minutes, but it should be run for every D/t')
tic;
%SymbolicVelocityDueToTube(D_over_T);
toc
disp('++Generate ODE source term::DONE')


%********************************************************************
%********************************************************************
%********************** VALIDATION AREA *****************************
%************* compare the CSL with the Closed form *****************
if ( true)
    figure(2105)
    zz = linspace(-20, 20, 400);
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




[RR, ZZ] = meshgrid([linspace(1e-6, 1, 60), linspace(1.05, 4, 60)],  ...
    linspace(6,-3,150) );
% allocatate memory
rr = 0*RR; zz = rr; disp_radial = rr; disp_vertical = rr;
eR = rr; eZ = rr; eTheta = rr; eRZ = rr;
vr = rr; vz = rr;

InitialZOfTube = -25;
FinalZOfTube = 3;




for jj = 1:size(RR,2)
    tic
    for ii = 1:size(RR,1)
        R = RR(ii,jj); Z = ZZ(ii,jj);
        
        [r, z, epsilon] = IntegrateDisplacements( InitialZOfTube, ...
            FinalZOfTube, R, Z);
        
        rr(ii,jj) = r; zz(ii,jj) = z;
        eR(ii,jj) = epsilon(1); eZ(ii,jj) = epsilon(2);
        eTheta(ii,jj) = epsilon(1); eRZ(ii,jj) = epsilon(4);
        [vr(ii,jj), vz(ii,jj)] = EvaluateVelocity(r, z, FinalZOfTube);
    end
    
    pause(0.0001)
    disp_radial = rr - RR;
    disp_vertical = zz - ZZ;
    
    save('StrainPath.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
        'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
        'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')
    if ( jj > 3)
        PostProcessResults('StrainPath', 'StrainPath.mat');
        disp(['PostProcess ', num2str(jj), ' of ', num2str( size(RR,2) )] )
        toc
    end
end


save('StrainPath.mat', 'rr','RR','zz','ZZ', 'disp_radial', ...
    'disp_vertical', 'ii', 'jj', 'vr', 'vz', ...
    'eR', 'eZ', 'eTheta', 'eRZ', 'D_over_T')



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


function [vr, vz] = EvaluateVelocity(r, z, FinalZ)

[a,~] = SourceTermStrainPath(FinalZ, r, z);
vr = a(1);
vz = a(2);

function [dxdt] = AuxiliarFunction(t, x)
r = x(1);
z = x(2);
[a,b] = SourceTermStrainPath(t,r,z);

dxdt = [a';b'];
