function [p] = pat_phantom(h, len_nodes)

ellipse = modified_shepp_logan;

p = zeros(len_nodes, 1);

x = 4 * ( h.fem.Promoted.nodes(1, :)');
y = 4 * ( h.fem.Promoted.nodes(2, :)');

for k = 1 : size(ellipse, 1)
   asq = ellipse(k,2)^2;       % a^2
   bsq = ellipse(k,3)^2;       % b^2
   phi = ellipse(k,6)*pi/180;  % rotation angle in radians
   x0 = ellipse(k,4);          % x offset
   y0 = ellipse(k,5);          % y offset
   A = ellipse(k,1);           % Amplitude change for this ellipse
   cosp = cos(phi);
   sinp = sin(phi);
   
   idx = find((((x - x0).*cosp + (y - y0).*sinp).^2)./asq + ...
       (((y - y0).*cosp - (x - x0).*sinp).^2)./bsq <= 1);
   p(idx) = p(idx) + A;
    
end


   

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Default head phantoms:   
%

function shep=shepp_logan
%
%  This is the default head phantom, taken from AK Jain, 439.
%
%         A    a     b    x0    y0    phi
%        ---------------------------------
shep = [  1   .69   .92    0     0     0   
        -.98 .6624 .8740   0  -.0184   0
        -.02 .1100 .3100  .22    0    -18
        -.02 .1600 .4100 -.22    0     18
         .01 .2100 .2500   0    .35    0
         .01 .0460 .0460   0    .1     0
         .01 .0460 .0460   0   -.1     0
         .01 .0460 .0230 -.08  -.605   0 
         .01 .0230 .0230   0   -.606   0
         .01 .0230 .0460  .06  -.605   0   ];
      
      
function toft=modified_shepp_logan
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A    a     b    x0    y0    phi
%        ---------------------------------
toft = [  1   .69   .92    0     0     0   
        -.8  .6624 .8740   0  -.0184   0
        -.2  .1100 .3100  .22    0    -18
        -.2  .1600 .4100 -.22    0     18
         .1  .2100 .2500   0    .35    0
         .1  .0460 .0460   0    .1     0
         .1  .0460 .0460   0   -.1     0
         .1  .0460 .0230 -.08  -.605   0 
         .1  .0230 .0230   0   -.606   0
         .1  .0230 .0460  .06  -.605   0   ];




