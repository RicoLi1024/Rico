
[p0,dir0,nouse]=xlsread('hk00005.xls');
[N M]=size(p0);

for i=1:N
  P(i)=p0(i,4);  

end
i=0;

for i=5:N
    pa=0;
    for j=1:3
        
     pa=pa+P(i-j+1)/3;   
    end
    p1(i)=pa;
end
i=0;
j=0;


for i=1:500
    Pa(i,:)=p1(5*i:492+5*i);
end
i=0;


for i=1:500
    P1(i,:)=P(5*i:492+5*i);
end
i=0;


r(:,1)=zeros(500,1);
for j=1:492 
    r(:,j+1)=log(P1(:,j+1)./P1(:,j));
end
j=0;


lmd=0.95;
p=[10 0;0 10];
c=0.01;
x3=log(P1./Pa);
for i=1:500
for k=1:492
    aa(:,k,i)=[0;0];
end;
end

for i=1:500
for k=1:492
    if abs(r(i,k))<0.1
    
    y1=meb(x3(i,k),0,c,2*c);
    y2=meb(x3(i,k),c,2*c,3*c);
    y3=meb(x3(i,k),2*c,3*c,3*c);
    y4=meb(x3(i,k),-2*c,-c,0);
    y5=meb(x3(i,k),-3*c,-2*c,-c);
    y6=meb(x3(i,k),-3*c,-3*c,-2*c);
    y7=meb(x3(i,k),-c,0,c);
    y=y1+y2+y3+y7;
    ed1=0;
    if y~=0
        ed1=(0.1*y1+0.2*y2+0.4*y3)/y;
    end;
    y=y4+y5+y6+y7;
    ed2=0;
    if y~=0
        ed2=(0.1*y4+0.2*y5+0.4*y6)/y;
    end;
    x=[ed1;ed2];
    
%    x=[abs(x13)+x13;abs(x13)-x13];
    %此部分为aa(:,:,k)的参数估计方法，策略确实有效，为了保障个人策略专利，此处设置为黑盒子，待真正入职后分享
    
    end;
    if abs(r(i,k))>=0.1
        aa(:,k,i)=aa(:,k-1,i);
    end;
end;

end
i=0;k=0;

for i=1:500
for k=3:492
alt=0;

for j=1:3
    
   alt=alt+aa(:,k-j+1,i)/3;
end
aa1(:,k,i)=alt;
end
end
%follow the big buyer
ho=0;nb=0;ns=0;
for i=1:500
    for k=1:492
        if aa1(2,k,i)>0 & aa1(1,k,i)>=0 & ho==0
            nb=nb+1;
            buy(i,nb)=k+1;
            ho=1;
        end
        if aa1(2,k,i)<=0 & ho==1
            ns=ns+1;
            sel(i,ns)=k+1;
            ho=0;
          
        end
    end 
         if nb>ns
           sel(i,nb)=493;
           ho=0;
         end
         nb=0;
         ns=0;
end
i=0;
j=0;
[Q,W]=size(buy);
rall(:,1)=ones(Q,1);

for i=1:Q
    
   
for j=1:W
if buy(i,j)~=0
rall(i,:)=rall(i,:)*(P1(i,sel(i,j))/P1(i,buy(i,j)));
end

end

end

arall=mean(rall);
sd1=std(rall);
sr1=(arall-1)/sd1;
i=0;j=0;


Mrall(:,1)=ones(Q,1);
for i=1:Q
    
   
for j=1:W
if buy(i,j)~=0
Mrall(i,j+1)=Mrall(i,j)*(P1(i,sel(i,j))/P1(i,buy(i,j)));
end

end

end

for i=1:Q
    for j=1:W
        if Mrall(i,j)~=0
        C=max(Mrall(i,1:j));
        if C==Mrall(i,j)
            RetraceRatio(i,j)=0;
        else
             RetraceRatio(i,j)=(Mrall(i,j)-C)/C;
            
        end
        end
    end
    
end
MaxRetraceRatio=-min( RetraceRatio,[],2);
AMaxRetraceRatio=mean(MaxRetraceRatio);
i=0;j=0;

%buy and hold
for i=1:500
        rall2(i,:)=P1(i,493)/P1(i,1);
      
end

arall2=mean(rall2);
sd2=std(rall2);
sr2=(arall2-1)/sd2;
AMaxRetraceRatio1=1-min(rall2);
% ride the mood
i=0;k=0;

for i=1:500
for k=5:492
alt=0;

for j=1:5
    
   alt=alt+aa(:,k-j+1,i)/5;
end
aa2(:,k,i)=alt;
end
end

ho=0;nb=0;ns=0;
for i=1:500
    for k=1:492
        if aa2(2,k,i)+aa2(1,k,i)>0 &  ho==0
            nb=nb+1;
            buy2(i,nb)=k+1;
            ho=1;
        end
        if aa2(2,k,i)+aa2(1,k,i)<0 & ho==1
            ns=ns+1;
            sel2(i,ns)=k+1;
            ho=0;
          
        end
    end 
         if nb>ns
           sel2(i,nb)=493;
           ho=0;
         end
         nb=0;
         ns=0;
end
i=0;
j=0;
[Q2,W2]=size(buy2);
rall3(:,1)=ones(Q2,1);

for i=1:Q2
    
   
for j=1:W2
if buy2(i,j)~=0
rall3(i,:)=rall3(i,:)*(P1(i,sel2(i,j))/P1(i,buy2(i,j)));
end

end

end

arall3=mean(rall3);
sd3=std(rall3);
sr3=(arall3-1)/sd3;
i=0;j=0;

Mrall3(:,1)=ones(Q2,1);
for i=1:Q2
    
   
for j=1:W2
if buy2(i,j)~=0
Mrall3(i,j+1)=Mrall3(i,j)*(P1(i,sel2(i,j))/P1(i,buy2(i,j)));
end

end

end

for i=1:Q2
    for j=1:W2
        if Mrall3(i,j)~=0
        C=max(Mrall3(i,1:j));
        if C==Mrall3(i,j)
            RetraceRatio3(i,j)=0;
        else
             RetraceRatio3(i,j)=(Mrall3(i,j)-C)/C;
            
        end
        end
    end
    
end
MaxRetraceRatio3=-min( RetraceRatio3,[],2);
AMaxRetraceRatio3=mean(MaxRetraceRatio3);
i=0;j=0;
%follow trend
for i=60:N
    pa=0;
    for j=1:60
        
     pa=pa+P(i-j+1)/60;   
    end
    p60(i)=pa;
end
i=0;
j=0;
for i=5:N
    pa=0;
    for j=1:5
        
     pa=pa+P(i-j+1)/5;   
    end
    p5(i)=pa;
end
i=0;
j=0;
for i=1:500
    Pa60(i,:)=p60(5*i:492+5*i);
end
i=0;
for i=1:500
    Pa5(i,:)=p5(5*i:492+5*i);
end

ho=0;nb=0;ns=0;
for i=1:500
    for k=1:492
        if Pa5(i,k)>Pa60(i,k) & Pa60(i,k)~=0 & ho==0
            nb=nb+1;
            buy3(i,nb)=k+1;
            ho=1;
        end
        if Pa5(i,k) < Pa60(i,k) & ho==1
            ns=ns+1;
            sel3(i,ns)=k+1;
            ho=0;
          
        end
    end 
         if nb>ns
           sel3(i,nb)=493;
           ho=0;
         end
         nb=0;
         ns=0;
end
i=0;
j=0;
[Q3,W3]=size(buy3);
rall4(:,1)=ones(Q3,1);

for i=1:Q3
    
   
for j=1:W3
if buy3(i,j)~=0
rall4(i,:)=rall4(i,:)*(P1(i,sel3(i,j))/P1(i,buy3(i,j)));
end

end

end

arall4=mean(rall4);
sd4=std(rall4);
sr4=(arall4-1)/sd4;
i=0;j=0;

Mrall4(:,1)=ones(Q3,1);
for i=1:Q3
    
   
for j=1:W3
if buy3(i,j)~=0
Mrall4(i,j+1)=Mrall4(i,j)*(P1(i,sel3(i,j))/P1(i,buy3(i,j)));
end

end

end

for i=1:Q3
    for j=1:W3
        if Mrall4(i,j)~=0
        C=max(Mrall4(i,1:j));
        if C==Mrall4(i,j)
            RetraceRatio4(i,j)=0;
        else
             RetraceRatio4(i,j)=(Mrall4(i,j)-C)/C;
            
        end
        end
    end
    
end
MaxRetraceRatio4=-min( RetraceRatio4,[],2);
AMaxRetraceRatio4=mean(MaxRetraceRatio4);
i=0;j=0;
t=linspace(1,500,500);
rmin1=min([rall rall2 rall3 rall4]);
rmin=min(rmin1);
rmax1=max([rall rall2 rall3 rall4]);
rmax=max(rmax1);
plot(t,rall,'Color','r','LineWidth',1.5);
hold on
plot(t,rall2,'Color','b','LineWidth',1.5);
hold on
plot(t,rall3,'Color','g','LineWidth',1.5);
hold on
plot(t,rall4,'Color','y','LineWidth',1.5);
xlabel('502个测试区间','FontSize',15,'FontWeight','bold');
ylabel('年化收益率（%）','FontSize',15,'FontWeight','bold');
legend('跟踪大买家策略','买入和持有策略','骑风而行策略','趋势跟踪策略','Location','best');
title('HK00005','Fontsize',12); 
ylim([rmin rmax]);
xlim([0 500]);
