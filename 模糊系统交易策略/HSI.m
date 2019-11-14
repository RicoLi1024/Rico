[p0,dir0,nouse]=xlsread('HSI.xls');
[N M]=size(p0);
for i=1:N
  P(i)=p0(i,4);  

end

t=linspace(1,N,N);

plot(t,P,'Color','b','LineWidth',1.5);

xlim([1 N]);