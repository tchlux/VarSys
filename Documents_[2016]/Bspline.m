bspl[j_Integer,1,t_List,x_Real] := If[t[[j]] <= x
< t[[j+1]], 1.0,0.0]

bspl[i_Integer,k_Integer,t_List,x_Real] := If[
t[[i+k-1]] == t[[i]], 0.0, (x-t[[i]])/(t[[i+k-1]] - t[[i]])*
bspl[i,k-1,t,x]] + If[t[[i+k]] == t[[i+1]], 0.0, (t[[i+k]]-x)
/(t[[i+k]]-t[[i+1]])*bspl[i+1,k-1,t,x]]

