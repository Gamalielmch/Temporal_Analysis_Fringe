% función peaks de 50x50
% Imagen de 50x50 y N=50, por lo tanto (t-1) t=1,2,3,...,N
% cycles ranging from 6 to 12
% clear all
% close all
mage=0:.005:0.035;
nimar=80;
for lo=1:length(nimar)

ni=nimar(lo);
f=peaks(100);
delta_phi_min=0.5010*pi;%pi/6; 0.5010*pi; (1/4)*pi;
delta_phi_max=0.8970*pi;%(1/2)*pi;%pi/2;pi*(2/3);
zmax=10;
zmin=-10;
dx=size(f,1);
miu=((delta_phi_max-delta_phi_min))/(zmax-zmin+dx);
x0=((delta_phi_min)/miu)-zmin;
a=128;
b=128;
ciclos = 6;
x=1:size(f,1);
X=size(f,1);
N=length(ciclos);
p0=X/ciclos(1);
fr=zeros(N,1);
cosenod=zeros(size(f));
ep1=(2*pi)*(1/p0);
% figure, colormap gray

for i=1:ni
    for j=1:size(f,2)
     cosenod(j,:,i)=a+b*cos(ep1*(f(j,:)+x0+x)+ miu*(f(j,:)+x0+x)*(i-1)+    (-0.03 + (0.03+0.03)*rand(1,100)));%+  (-0.8 + (0.8+0.8)*rand(1,100));
    end
%       imagesc(coseno(1,:,i))
%       colormap(gray)
end


% r = -0.8 + (0.8+0.8)*rand(100,1);


W=wavMadre(21,6);
freco=zeros(size(f));
% freco2=freco;
Ns =21;
Nt=13;
t=6;
lon=ni-2*t;
rid=zeros(lon,1);
Tm=1:ni;
Tm=Tm(t+1:ni-t);
m=length(Tm);
na=(m*sum(Tm.^2)-sum(Tm)^2)';
na2=(m*Tm)';
na3=sum(Tm)';
M=size(freco,1);
N=size(freco,2);
WT=zeros(Ns,ni-2*t);
ve=0:21:21*(ni-Nt);
[x,~]=meshgrid(1:13,1:ni-Nt+1);
[j2,~]=meshgrid(0:size(freco,2)-1,0:size(freco,1)-1);
for i=2:(ni-Nt+1)
x(i,:)=x(i,:)+i-1;
end
%h = waitbar(0,'Please wait...');
le=size(freco,1)*size(freco,2);
tic
parfor i=1:size(freco,1)*size(freco,2)
        %waitbar(i/le,h)
        te=cosenod(i+(M*N*(0:ni-1)))';
        te2=te(x);
        WT=W*te2';
        [mx,cm]=max(abs(WT));
        rid=WT(cm+ve);
        %Calcula la frecuencia de la señal
        ph=unwrap(angle(rid));
        freco(i)=(ph*na2-na3*sum(ph))/na;
%         freco2(i)=d1sincrona(te',1:length(te));
end
freco=(freco/miu)-j2;
% freco2=(freco2/miu)-j2;
tie(lo) = toc;
of=min(min(freco))-min(min(f));
freco=freco-of;
% of=min(min(freco2))-min(min(f));
% freco2=freco2-of;
% figure(1)
% mesh(freco);
% [er1(lo), er2(:,:,nii)]=rmserror(f,freco);
[er1(lo), ~]=nrmserror(f,freco);
% [er2(lo),~]=nrmserror(f,freco2);
% figure
% mesh(er2)
end