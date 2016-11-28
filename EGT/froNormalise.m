function F = normaliseF(Fin)
% Normalises F so 3,3 = 1 but keeps its vec or matrix shape%newLength = size(Fin,1)*size(Fin,2);%vecF = reshape(Fin,1,newLength);
s = sign(Fin(size(Fin,1),size(Fin,2)));
if s == 0    s = 1;endF = Fin/(s*norm(Fin,'fro'));

    