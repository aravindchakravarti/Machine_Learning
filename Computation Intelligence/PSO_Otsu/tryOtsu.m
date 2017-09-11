% --------- PARTICLE SWARM OPTIMIZATION FOR IMAGE SEGMENTATION ------------
%
% File Name     : tryOtsu.m
%
% Description   : This function is implementation of OTSU function. This
%                 function takes inputs and returns the fitness value
%                 achieved for those inputs.
%                 Higher the fitness, maximum is the 'Betwen Class
%                 Variance' of image.
%
% Reference     : A Threshold Selection Method from Gray-Level Histograms 
%                 IEEE Transactions on Systems, Man, and Cybernetics, 
%                 Vol.9, No. 1, 1979, pp. 62-66.
%                 Author: Noboyuki Otsu
%
% Function      : image_hist_p - Normalized image histogram
% Inputs          no_level     - How many levels thresholding required
%                 k            - Particles representing possible thresholds
%                 mau_t        - Image mean
%
% Function      : fitness      - Fitness value acheived
% Outputs
%
% File History  : | Ver. |  Date     |              Changes               |
%                 | 1.0  | 03/AUG/17 | Base version of Otsu fucntion      |
%                 |      |           |                                    |
%                 | 1.1  | 03/SEP/17 | Comments are added                 |
%                 |      |           |                                    |
%                 |      |           |                                    |
%                 |      |           |                                    |
% -------------------------------------------------------------------------
function fitness = tryOtsu(image_hist_p, no_level, k, mau_t)
    
    % Zero'th and First'th order moment
    omega = zeros (1, no_level + 1);
    mau   = zeros (1, no_level + 1);
    
    % For number levels
    for level_indx = 1:(no_level+1)
        if level_indx == 1
            start_pos = 1;
            end_pos   = k(level_indx);
        elseif ((level_indx > 1) && (level_indx < (no_level+1)))
            start_pos = k(level_indx-1)+1;
            end_pos   = k(level_indx);
        else
            start_pos = k(level_indx-1)+1;
            end_pos   = 256;    
        end
        
        % Calculation of zeroth and first order moment
        omega(level_indx) = sum(image_hist_p(start_pos:end_pos));
            
        for i = start_pos: end_pos
            mau(level_indx) = mau(level_indx) + (i*image_hist_p(i)...
                                                    /omega(level_indx));
        end
    end
         
    sigma_b_sqr = 0;
    
    % This is the calculation of between class variance
    for i=1:(no_level+1)
        sigma_b_sqr = sigma_b_sqr + (omega(i).*((mau(i)-mau_t).^2));
    end
    
    fitness = sigma_b_sqr;
   
end
