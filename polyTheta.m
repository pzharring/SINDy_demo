function Theta = polyTheta(x,order)

len = size(x,1);
dims = size(x,2);
colnum = 1;

Theta(:,colnum) = ones(len,1);
colnum = colnum + 1;

if order >= 1
    for i=1:dims
        Theta(:,colnum) = x(:,i);
        colnum = colnum + 1;
    end
end

   
if order >= 2
    for i = 1:dims
        for j = i:dims
            Theta(:,colnum) = x(:,i).*x(:,j);
            colnum = colnum + 1;
        end
    end
end

if order >= 3
    for i = 1:dims
        for j = i:dims
            for k = j:dims
                Theta(:,colnum) = x(:,i).*x(:,j).*x(:,k);
                colnum = colnum + 1;
            end
        end
    end
end

if order==4
    for i = 1:dims
        for j = i:dims
            for k = j:dims
                for l = k:dims
                    Theta(:,colnum) = x(:,i).*x(:,j).*x(:,k).*x(:,l);
                    colnum = colnum + 1;
                end
            end
        end
    end
end
