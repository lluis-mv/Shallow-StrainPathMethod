%=========================================================================
%  Symbolic Matlab to compute the velocities and strains Strains due to
%  the penetration of a tube using the Strain Path Method (Baligh, 1985)
%
%
%  Following the approach of Sagaseta et al (1991)
%
%  Barcelona, 13 February 2019
%=========================================================================

function [velocity, v1, v2, v3]=SourceTermShallowStrainPath(h, r, z)

% 1- Define the geometry of the tube
global D_over_T;
if ( size(D_over_T, 1) == 0)
    D_over_T = 40;
end

if ( nargin == 0)
    h = 0.1;
    r = 2;
    z = 3;
end

R = 1;
t = 2*R/D_over_T;
U = 1;


% Sagaseta et al (1991) define differently the geometry of the problem
%    see Figure 3 of the referred work
w = t/2;
R = R-w;

%**************************************************************************
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
vr1 = U * w  * R / pi * ( A1rK*K + A1rE * E );
vz1 = U * w * R / pi * ( A1zE * E);

%**************************************************************************
% Appendix A, Equation (29) of Sagaseta et al (1991)
r21_pow2 = ( (r-R)^2 + (z+h)^2);
r21 = sqrt(r21_pow2);
r22_pow2 = ( (r+R)^2 + (z+h)^2);
r22 = sqrt(r22_pow2);

k2_pow2 =  1 - r21_pow2/ r22_pow2 ;


%[K,E] = ellipke(k1_pow2);
[K,E ] = ApproximateKandE(k2_pow2);
A2rK = 1/r22 / r;
A2rE = -1/r22*(1/r - 2 * ( r-R)/r21_pow2 );
A2zE = 1/r22 * 2 * (z+h)/r21_pow2;

% Velocities due to a tube
vr2 = -U * w  * R / pi * ( A2rK*K + A2rE * E );
vz2 = -U * w * R / pi * ( A2zE * E);

%**************************************************************************
% Appendix A, Equation (33) of Sagaseta et al (1991)
[vr3, vz3] = EvaluateVelocity3( h, r, z, R, w);
% vr3 = 0; vz3 = 0;





velocity = [vr1+vr2+vr3, vz1+vz2+vz3]';
if ( nargout > 1)
    v1 = [vr1, vz1]';
    v2 = [vr2, vz2]';
    v3 = [vr3, vz3]';
end


function [vr3, vz3] = EvaluateVelocity3( h, r, z, R, w)


fun = @(x) SourceTermIntegral(x, h, r, z, R, w, 'R');
vr3 = integral( fun, 0, Inf, 'AbsTol', 6e-8);

fun = @(x) SourceTermIntegral(x, h, r, z, R, w, 'Z');
vz3 = integral( fun, 0, Inf, 'AbsTol', 6e-8);

function [ss] = SourceTermIntegral( a, h, r, z, R, w, COMPONENT)


rh1_pow2 = (a-R).^2 + h^2;
rh1 = sqrt(rh1_pow2);
rh2_pow2 = (a+R).^2 + h^2;
rh2 = sqrt(rh2_pow2);

kh_pow2 = 1 - rh1_pow2./rh2_pow2;
[K,E ] = ApproximateKandE(kh_pow2);

AhK = - h./(rh1_pow2.*rh2) .* ( 1-kh_pow2)./kh_pow2 .* ...
        (  (a-R)./rh1_pow2 - (a+R)./rh2_pow2);

AhE = 2*(a-R)./rh1_pow2   + (a+R)./rh2_pow2 + ( 1-kh_pow2)./kh_pow2 .* ...
    ( (a-R)./rh1_pow2 - (a+R)./rh2_pow2);

AhE = AhE * h ./ rh1_pow2./ rh2;

FgammaST = 8 * 1 * w * R / pi * ( AhK .* K + AhE .* E);

ra1_pow2 = (a-r).^2 + z^2;
ra1 = sqrt(ra1_pow2);
ra2_pow2 = (a+r).^2 + z^2;
ra2 = sqrt(ra2_pow2);

k_pow2 = 1 - ra1_pow2./ra2_pow2;
[K,E ] = ApproximateKandE(k_pow2);

ArK = ra2/r.*( 1 - (2*a*r- z^2)./ra2_pow2);
ArE = -ra2./r.*( 1 + z^2/2*( 1./ra1_pow2 + 1./ra2_pow2) );
AzK = -z./ra2;
AzE = z./ra2 .* ( 1 - 2*a.*(a-r)./ra1_pow2);

Ir = 1/2/pi * ( ArK.*K + ArE.*E);
Iz = 1/2/pi * ( AzK.*K + AzE.*E);

sr = Ir.*FgammaST;
sz = Iz.* FgammaST;

ss = sz;
if ( COMPONENT == 'R')
    ss = sr;
end

function b = myFactorial(n)
b = prod(1:n);


function [K, E] = ApproximateKandE(x2, nTerms)
% http://www.exstrom.com/math/elliptic/ellipint.html
if (nargin == 1)
    nTerms = 70;
end

K = 0;
for n = 0:nTerms
    number = (myFactorial(2*n)/myFactorial(n)/myFactorial(2*n-n))^2 / (4^(2*n));
    K = K + number * x2.^n;
end
E = 1;
for n = 1:nTerms
    number = (myFactorial(2*n)/myFactorial(n)/myFactorial(2*n-n))^2 / (4^(2*n)) / ( 2*n-1);
    E = E - number * x2.^n ;
end
K = K*pi/2;
E = E*pi/2;