% W = Familia de wavelets
% Ne = Número de escalas
% tt = Define el tamaño de la wavelet: 2*tt+1
% Se recomienda Ne=21 y tt=6 para experimentos de perfilometría con 80
% imágenes

function W=wavMadre(Ne,tt)

N=2*tt+1;
W=zeros(Ne,N);
x=[-tt:tt];
frec=linspace(1/7,1/2,Ne);
sig=linspace(7,2,Ne);

for k=1:Ne
    W(k,:)=exp(-i*2*pi*frec(k)*x).*exp(-0.5*(x/sig(k)).^2);
end    
    


