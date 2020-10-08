% función peaks de 50x50
% Imagen de 50x50 y N=50, por lo tanto (t-1) t=1,2,3,...,N
% cycles ranging from 6 to 12
% clear all
% close all
nii=1;
nima=80;
f=abs(peaks(100));
zmax=10;
zmin=-10;
dx=size(f,1);
xt=1:size(f,1);
p0=size(f,1)/6;
cosenod=zeros(size(f));
ep1=(2*pi)*(1/p0);
minimo=300;
tic
delta_phi_min=pi*(0.5000:.0001:0.5020);%pi/6; %0.94*pi;
delta_phi_max=pi*(0.8900:.001: 0.90);%pi/2; %pi*0.31;
W=wavMadre(21,6);
freco=zeros(size(f));
Ni=80;
Ns =21;
Nt=13;
t=6;
lon=Ni-2*t;
rid=zeros(lon,1);
Tm=7:74;
m=length(Tm);
na=(m*sum(Tm.^2)-sum(Tm)^2)';
na2=(m*Tm)';
na3=sum(Tm)';
M=size(freco,1);
N=size(freco,2);
WT=zeros(Ns,Ni-2*t);
ve=0:21:21*67;
[x,~]=meshgrid(1:13,1:68);
[j2,~]=meshgrid(0:size(freco,2)-1,0:size(freco,1)-1);
for i=2:68
    x(i,:)=x(i,:)+i-1;
end
er1=zeros(length(delta_phi_min),length(delta_phi_max));
er3=er1;


for nn=1:length(delta_phi_min)
    for mm=1:length(delta_phi_max)
        if delta_phi_max(mm)>delta_phi_min(nn)
            miu=((delta_phi_max(mm)-delta_phi_min(nn)))/(zmax-zmin+dx);
            x0=((delta_phi_min(nn))/miu)-zmin;
            try
                
                for i=1:nima
                    for j=1:size(f,2)
                        cosenod(j,:,i)=128+128*cos(ep1*(f(j,:)+x0+xt)+ miu*(f(j,:)+x0+xt)*(i-1));%+  (-0.8 + (0.8+0.8)*rand(1,100));
                    end
                    %       imagesc(coseno(1,:,i))
                    %       colormap(gray)
                end
                
                parfor i=1:size(freco,1)*size(freco,2)
                    %waitbar(i/le,h)
                    te=cosenod(i+(M*N*(0:79)))';
                    te2=te(x);
                    WT=W*te2';
                    [mx,cm]=max(abs(WT));
                    rid=WT(cm+ve);
                    %Calcula la frecuencia de la señal
                    ph=unwrap(angle(rid));
                    freco(i)=(ph*na2-na3*sum(ph))/na;
                end
                freco=(freco/miu)-j2;
                of=min(min(freco))-min(min(f));
                freco=freco-of;
                %         mesh(freco);
                [er1(nn,mm), er2]=rmserror(f,freco);
                [er3(nn,mm), er4]=nrmserror(f,freco);
                if er3(nn,mm)<minimo
                    minimo=er3(nn,mm);
                    pmin=delta_phi_min(nn);
                    pmax=delta_phi_max(mm);
                    buena=freco;
                end
            catch
                er1(nn,mm)=100; er3(nn,mm)=100;
            end
        end
    end
end

toc
