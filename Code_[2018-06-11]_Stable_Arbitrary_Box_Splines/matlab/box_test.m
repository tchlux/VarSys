%% k = 3 directions in R^s = R^2. Each multiplicity two.

echo off

tic
 ## X      = [ 1  0  1 ;
 ##            0  1  1 ];
 ## nu     = [2;2;2];
 X      = [ 1  0   1  1 ;
            0  1  -1  1 ];
 nu     = [2;2;2;2];
## [xx,yy] = meshgrid(((1:20)-2)/5,((1:20)-2)/5);
[xx,yy] = meshgrid((1:50),(1:50));
 p      = [xx(:) yy(:)];
disp(size(p))
toc

tic
 ## p      = p(1:5,:);
 ## p = [ 1.0001 1; ];
 b = box_eval(X,nu,p);
 ## disp('b:')
 ## disp(b)
 ## surf(reshape(b,20,20))
 ## pause
toc

## 58.433 seconds total
