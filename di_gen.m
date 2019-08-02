

function im_di = di_gen(im1, im2)
    im1 = srad( im1, 0.15, 1, 1, 1, 5 );
    im2 = srad( im2, 0.15, 1, 1, 1, 5 );
    im_di = abs(log( (im1+1) ./ (im2+1) ));
    im_di = srad( im_di, 0.15, 1, 1, 1, 5 );
end


function im = srad( im,delta_t,q0,rho,option,ITER ) 
    for n = 1:ITER 
        % padding elements avoid edges
        [m, n] = size(im);
        W = 1;
        X = zeros(m+2*W, n+2*W);
        X(W+1:m+W, W+1:n+W) = im;
        X(W+1:m+W, 1:W) = im(:, W:-1:1);
        X(W+1:m+W, n+W+1:n+2*W) = im(:, n:-1:n-(W-1));
        X(1:W, :) = X(2*W:-1:(W+1), :);
        X(m+(W+1):m+2*W, :) = X(m+W:-1:(m+1), :);
        im = X;
        % padding finished !!!!

        q0 = q0 * exp(-delta_t*rho); 

        gradRxI = imfilter(im,[0 0 0;0 -1 1;0 0 0],'symmetric'); 
        gradRyI = imfilter(im,[0 1 0;0 -1 0;0 0 0],'symmetric'); 
        gradLxI = imfilter(im,[0 0 0;-1 1 0;0 0 0],'symmetric'); 
        gradLyI = imfilter(im,[0 0 0;0 1 0;0 -1 0],'symmetric'); 
        q1 = sqrt(gradRxI.^2+gradRyI.^2+gradLxI.^2+gradLyI.^2)./( im+1e-10 ); 
        q2 = 4*del2(im)./(im+1e-10); 
        q = sqrt((1/2*(q1.^2)-1/16*(q2.^2))./((1+1/4*q2).^2+1e-10)); 

        if option == 1 
            c = 1./(1+((q.^2-q0^2)/(q0^2*(1+q0^2)))); 
        elseif option == 2 
            c = exp(-(q.^2-q0^2)/(q0^2*(1+q0^2))); 
        else 
            error('Unidentified option!(only 1 or 2 available)'); 
        end 

        d = imfilter(c,[0 0 0;0 0 1;0 0 0],'symmetric').*... 
            imfilter(im,[0 0 0;0 -1 1;0 0 0],'symmetric')+...         
            imfilter(c,[0 0 0;0 1 0;0 0 0],'symmetric').*... 
            imfilter(im,[0 0 0;1 -1 0;0 0 0],'symmetric')+... 
            imfilter(c,[0 1 0;0 0 0;0 0 0],'symmetric').*... 
            imfilter(im,[0 1 0;0 -1 0;0 0 0],'symmetric')+... 
            imfilter(c,[0 0 0;0 1 0;0 0 0],'symmetric').*... 
            imfilter(im,[0 0 0;0 -1 0;0 1 0],'symmetric'); 

        im = im+delta_t/4*d; 

        im = im(2:end-1,2:end-1); 
    end
end





