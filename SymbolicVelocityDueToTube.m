%=========================================================================
%  Symbolic Matlab to compute the velocities and strains Strains due to
%  the penetration of a tube using the Strain Path Method (Baligh, 1985)
%
%
%  Following the approach of Sagaseta et al (1991)
%
%  Barcelona, 13 February 2019
%=========================================================================

function []=SymbolicVelocityDueToTube()

% 1- Define the geometry of the tube
global D_over_T;
if ( size(D_over_T, 1) == 0)
    D_over_T = 40;
end

syms h positive;

R = 1;
t = 2*R/D_over_T;
U = 1;

% 2- Coordinates of an spatial point
syms r positive;
syms z real;

% Sagaseta et al (1991) define differently the geometry of the problem
%    see Figure 3 of the referred work
w = t/2;
R = R-w;

% Appendix A, Equation (25) of Sagaseta et al (1991)
r11_pow2 = ( (r-R)^2 + (z-h)^2);
r11 = sqrt(r11_pow2);
r12_pow2 = ( (r+R)^2 + (z-h)^2);
r12 = sqrt(r12_pow2);

k1_pow2 =  1 - r11_pow2/ r12_pow2 ;

% Matlab and Sagaseta et al (1991) use different difinitions of the
% elliptical integrals
%[K,E] = ellipke(k1_pow2);

[K,E ] = ApproximateKandE(k1_pow2);
A1rK = 1/r12 / r;
A1rE = -1/r12*(1/r - 2 * ( r-R)/r11_pow2 );
A1zE = 1/r12 * 2 * (z-h)/r11_pow2;

% Velocities due to a tube
vr = U * w  * R / pi * ( A1rK*K + A1rE * E );
vz = U * w * R / pi * ( A1zE * E);

% Temporal derivative of the strains due to a tube
epsiZdot = -diff( vz, z);
epsiRdot = -diff(vr, r);
epsiThetadot = -vr/r;
epsiRZdot = -1/2*( diff(vr, z) + diff(vz, r));


% Export the velocities and temporal derivative of strains to a file
%      to be used in the integration to obtain displacements and strains
velocity = [vr, vz];
epsiDot = [epsiRdot, epsiZdot, epsiThetadot, epsiRZdot];

matlabFunction(velocity, epsiDot, 'File', 'SourceTermStrainPath', 'Vars', [h, r, z]);

return;
%********************************************************************
%********************************************************************
%********** debug to check that the derivatives are correct *********
if ( false)
    
    nu1 = k1_pow2;
    dnu1dr = diff(nu1,r);
    dnu1dz = diff(nu1,z);
    
    eTheta = -U*w*R/pi/r * ( A1rK*K + A1rE*E);
    eZZ = -U*w*R/pi*( - A1zE/(2*nu1)*dnu1dz*K + (diff(A1zE,z) + A1zE/(2*nu1)*dnu1dz )* E );
    ezr = -U*w*R/pi*( - A1zE/(2*nu1)*dnu1dr*K + (diff(A1zE,r) + A1zE/(2*nu1)*dnu1dr )* E );
    
    ErrorInEpsiT = simplify(epsiThetadot - eTheta)
    ErrorInEpsiZ = simplify(epsiZdot - eZZ)
    
    ErrorInEpsiZR = simplify(epsiRZdot - ezr)
    ErrorInEpsiVol = simplify(epsiZdot + epsiRdot + epsiThetadot)
    
end
%********************************************************************
%********************************************************************


%********************************************************************
%********************************************************************
%********** check that the apKnumproximation of elliptical **************
% ***************  integrals is correct *****************************
if ( true)
    syms k_pow2 positive
    
%     [Kanal,Eanal] = ellipke(k_pow2);
%     
%     for nTerms = [4:4:32]
%         [Knum, Enum] = ApproximateKandE(k_pow2, nTerms);
%         
%         errorK = sqrt( (Kanal-Knum)^2);
%         errorE = sqrt( (Eanal-Enum)^2);
%         
%         x = linspace(1e-6,1, 100);
%         err1 = x;
%         err2 = x;
%         
%         for i = 1:length(x)
%             err1(i) = eval(subs( errorK, k_pow2, x(i)));
%             err2(i) = eval(subs( errorE, k_pow2, x(i)));
%         end
%         figure(121)
%         semilogy( x, err1)
%         hold on
%         semilogy(x, err2)
%         hold off;
%         pause(1)
%     end
%     
    
    % Second part
    nTerms = 2:2:100;
    [Kanal,Eanal] = ellipke(k_pow2);
    x = 0.98;
    err1 = 0*nTerms;
    err2 = 0*nTerms;
    
    for i = 1:length(nTerms)
        nn = nTerms(i);
        
        [Knum, Enum] = ApproximateKandE(k_pow2, nn);
        
        errorK = sqrt( (Kanal-Knum)^2);
        errorE = sqrt( (Eanal-Enum)^2);
        
        err1(i) = eval(subs( errorK, k_pow2, x));
        err2(i) = eval(subs( errorE, k_pow2, x));
        
        figure(122)
        loglog( nTerms, err1)
        hold on
        semilogy(nTerms, err2)
        hold off;
        drawnow
    end
    
    
end


function [K, E] = ApproximateKandE(x2, nTerms)
% http://www.exstrom.com/math/elliptic/ellipint.html
if (nargin == 1)
    nTerms = 50;
end


K = 0;
for n = 0:nTerms
    n1 = sym(n);
    number = (factorial(2*n1)/factorial(n1)/factorial(2*n1-n1))^2 / (4^(2*n1));
    K = K + number* x2^n; 
end
E = 1;
for n = 1:nTerms
    n1 = sym(n);
    number = (factorial(2*n1)/factorial(n1)/factorial(2*n1-n1))^2 / (4^(2*n1)) / ( 2*n1-1);
    E = E - number * x2^n ;
end

K = K*pi/2;
E = E*pi/2;