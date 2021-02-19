function simulacion2V
w1=(2*pi)/5;
w2=(2*pi)/2.5;
N=50;
sz=100;
theta=45*(pi/180);
z=peaks(sz);
ws=(w2-w1)/sz;
wt=(meshgrid(0:sz-1)*cos(theta)+z*sin(theta))*ws+w1;
I=zeros(sz,sz,N);
wtr=zeros(sz,sz);
for t=1:N
        I(:,:,t)=1+cos(wt*t);
end
Na=50;
es=linspace(2,7,Na);
es2=repmat(es',1,N);
sg=1.75;
W=MoWav2(es,N,sg);
ran=round(N/2)-ceil(N*.1):round(N/2)+ceil(N*.1);
tic
for i=1:sz*sz
   te=I(i+(sz*sz*(0:N-1)))';
    te=repmat(reshape(te,1,N),Na,1,N);
    w=squeeze(sum(te.*W,2));
    Sap=(abs(w)./es2).^4;
    ar=sum((Sap.*es2),1)./sum(Sap,1);
    k=mean(ar(ran));
    wtr(i)=(2*pi)/k;
end
toc
figure, mesh(wtr-wt)
Nrmse=sqrt(abs(wt-wtr).^2./abs(wt).^2);
Ew=wt-wtr;
sqrt(sum(Ew(:).^2)/sum(wt(:).^2))
mean(Nrmse(:))
wtr=(((wtr-w1)/ws)-meshgrid(0:sz-1)*cos(theta))/sin(theta);
figure, mesh(wtr-z)
figure,mesh(wtr)
end



function W=MoWav2(es,N,sg)

Ne=length(es);
W=zeros(Ne,N,N);
x=(0:N-1)';
for a=1:Ne
    for b=0:N-1
        W(a,:,b+1)=exp(-1i*2*pi*(x-b)/es(a)).*exp(-((x-b)/es(a)).^2/(2*sg));
    end
end
end















