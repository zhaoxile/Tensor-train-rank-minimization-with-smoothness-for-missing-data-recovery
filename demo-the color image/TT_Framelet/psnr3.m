function p = psnr3( X, Y )         % pixel values of images are normalized to [0,1].

% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       psnr(X,Y) = 10*log10( max(max(X),max(Y))^2 / |X-Y|^2 ).
%
% p = psnr( X, Y );

sizeX = size( X );                 

if numel(sizeX) == 2               % compute the psnr value of the gray image
    
    d  = mean( mean( ( X(:)-Y(:) ).^2 ) );
    m1 = max( abs( X(:) ) );
    m2 = max( abs( Y(:) ) );
    m  = max( m1, 2 );

    p = 10*log10( m^2/d ); 
    
elseif numel(sizeX) == 3          % compute the average psnr value of the 3-order image
    
    psnr_vector = zeros(1,sizeX(end));
    
    for i = 1:sizeX(end)
        
        XX = X( :, :, i );
        YY = Y( :, :, i );
        d  = mean( mean( ( XX(:)-YY(:)).^2 ) );
        m1 = max( abs( XX(:) ) );
        m2 = max( abs( YY(:) ) );
        m  = max( m1, m2 );
        p  = 10*log10( m^2/d );
        psnr_vector(i) = p;
        
    end
    
    p = mean(psnr_vector);
    
else                             % output error when dimension is greater than 3.
    
    fprintf('Error');
    
end


