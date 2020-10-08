
w1=(2*pi)/5;
w2=(2*pi)/2.5;
N=30;
sx=512;
I2=zeros(sx,sx,N);
z=peaks(sx);
Fd=I2;
W=wavMadre(21,6);


for t=1:N;
ws=linspace(w1*t,w2*t,sx);
I=cos(ws);
I2(:,:,t)=repmat(I,sx,1);
for j=1:sx
Fd(j,:,t)=cos(ws+z(:,j)');
end
%imagesc(Fd(:,:,t))
%imagesc(I)
%title(num2str(t))
%colormap gray
%pause(0.3)
% print(['C:\Users\Gama\Documents\togif\', num2str(t),'.png'],'-dpng') %%% Change Path
end

[x,~]=meshgrid(1:13,1:68);
Ir=zeros(sx,sx);
Tm=7:74;
m=length(Tm);
na=(m*sum(Tm.^2)-sum(Tm)^2)';
na2=(m*Tm)';
na3=sum(Tm)';
M=size(z,1);
Nn=size(z,2);
ve=0:21:21*67;

for i=1:sx*sx
    %waitbar(i/le,h)
    te=Fd(i+(M*Nn*(0:N-1)))';
    te2=te(x);
    WT=W*te2';
    [mx,cm]=max(abs(WT));
    rid=WT(cm+ve);
    %Calcula la frecuencia de la señal
    ph=unwrap(angle(rid));
    Ir(i)=(ph*na2-na3*sum(ph))/na;
end
[j2,~]=meshgrid(0:size(z,2)-1,0:size(z,1)-1);
miu=1;
Ir=(Ir/miu)-j2;
figure, mesh(Ir)