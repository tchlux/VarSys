%% k = 3 directions in R^s = R^2. Each multiplicity two.

echo off

 ## X      = [ 1  0  1 ;
 ##            0  1  1 ];
 ## nu     = [2;2;2];
 X      = [ 1  0   1  1 ;
            0  1  -1  1 ];
 nu     = [1;1;1;1];

[xx,yy] = meshgrid(((1:20)-2)/5,((1:20)-2)/5);
 p      = [xx(:) yy(:)];
 ## p      = p(1:5,:);
 p = [ 1.0001 1; ];
 b = box_eval(X,nu,p);
 disp('b:')
 disp(b)
 ## surf(reshape(b,20,20))
 ## pause

