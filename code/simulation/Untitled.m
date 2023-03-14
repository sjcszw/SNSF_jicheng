
[U,S,V] = svd(tem3);
U*S*V'-tem3
%%
G = V*[inv(S(1:15,1:15)) rand(15,9);rand(30,24)]*U';
tem3*G*tem3-tem3
%%
G1 = V*[inv(S(1:15,1:15)) rand(15,9);zeros(30,24)]*U';
G2 = V*[inv(S(1:15,1:15)) rand(15,9);zeros(30,24)]*U';
tem3*G1 - tem3*G2
%%