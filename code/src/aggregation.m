function res=aggregation(P,s)
c=floor(s/2)+1;
n=size(P,1);
m=size(P,2);
isInside=@(x,y) (x>=1)&(x<=n)&(y>=1)&(y<=m);
mask=zeros(size(P,1),size(P,2));
res=zeros(size(mask));
    for i=1:size(P,1)
       for j=1:size(P,2)
           for k=1:size(P,3)
              for l=1:size(P,4)
                 if (isInside(i+k-c,j+l-c)) 
                     res(i+k-c,j+l-c)=res(i+k-c,j+l-c)+P(i,j,k,l); 
                     mask(i+k-c,j+l-c)=mask(i+k-c,j+l-c)+1;
                 end
              end
           end
       end
    end
res=res./mask;
end