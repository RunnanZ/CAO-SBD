%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aberrationDecomposition decomposes the aberration function into Zernike %
% coefficients                                                            %
%                                                                         %
% Inputs:                                                                 %
%        aberration   : aberration of the imaging system                  %
%        zernike_poly : Zernike polynomials                               %
% Outputs:                                                                %
%        zernike_coeff: Zernike coefficients for the aberration function  %
function zernike_coeff = aberrationDecomposition(aberration)
    global zernike_poly;
    zernike_coeff = zernike_poly'*aberration(:);%zernike_poly' Zt转置
end