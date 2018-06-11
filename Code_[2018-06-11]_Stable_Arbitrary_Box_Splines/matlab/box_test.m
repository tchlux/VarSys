%% k = 3 directions in R^s = R^2. Each multiplicity two.

echo off

 X      = [ 1  0  1 ;
            0  1  1 ];
 nu     = [2;2;2];
[xx,yy] = meshgrid(((1:20)-2)/5,((1:20)-2)/5);
 p      = [xx(:) yy(:)];

 ## p      = p(1:5,:);
 p = [ 2 2; ];
 b = box_eval(X,nu,p);

 ## surf(reshape(b,20,20))
 ## pause

