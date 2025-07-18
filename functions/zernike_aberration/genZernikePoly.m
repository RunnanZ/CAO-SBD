%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genZernikePoly generates a finite number of Zernike polynominals, which %                                                  %
% serves as the bases of pupil aberration function.                       %
%                                                                         %
% Inputs:                                                                 %
%        Fx           : spatial frequencies along x-axis on 2D grid       %
%        Fy           : spatial frequencies along y-axis on 2D grid       %
%        na           : the numerical aperture of the imaging system      %
%        lambda       : wavelength of the illumination                    %
%        highest_order: the highest order of Zernike bases to be evaluated%
% Outputs:                                                                %
%        zernike_poly : Zernike polynomials                               %
function zernike_poly = genZernikePoly(Fx, Fy, na, lambda, highest_order)
    [theta, rho] = cart2pol(Fx, Fy);
    rho          = rho./(na/lambda);
    zernike_poly = zeros(size(rho, 1)*size(rho,2), highest_order);
    pupil        = (Fx.^2+Fy.^2<=(na/lambda)^2);
    
    for zernike_index = 0:highest_order-1
        n       = ceil((-3+sqrt(9+8*zernike_index))/2);
        m       = 2*zernike_index - n*(n+2);
        zernike = zeros(size(rho));
        for k = 0:(n-abs(m))/2
            zernike = zernike + ((-1)^k*factorial(n-k))/...
                                (factorial(k)*factorial(0.5*(n+m)-k)*factorial(0.5*(n-m)-k))...
                                .*rho.^(n-2*k);
        end
        zernike = pupil.*zernike.*((m<0)*sin(abs(m)*theta) + (m>=0)*cos(abs(m)*theta));
        zernike_poly(:,zernike_index+1) = zernike(:);
    end
    
    % discard the first three Zernike polynomials
    zernike_poly = zernike_poly(:, 4:end);  
end