function xprime = genericpoly(tin,yin,Xid,nVars,inoned,polyord,...
    laurentord,absshifted,sineord,intimed,sintimed);
%Computes the coordinates in a generic polynomial equation.
%Repeats the poolData3 structure, but for single value inputs

m = size(yin,1); 
%yout is the alphabet matrix, here a row vector

ind = 1;

% poly order 0
if(polyord>=0)
    if (inoned==1)
         yout(ind) = 1;
         %yout(ind) = ones(m,1);
         ind = ind+1;
    end
end
% poly order 1
if(polyord>=1)
    for i=1:nVars
        yout(ind) = yin(i);
        ind = ind+1;
    end
end
% poly order 2
if(polyord>=2)
    for i=1:nVars
        for j=i:nVars
            yout(ind) = yin(i)*yin(j);
            ind = ind+1;         
        end
    end
end
% poly order 3
if(polyord>=3)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(ind) = yin(i)*yin(j)*yin(k);
                ind = ind+1;
            end
        end
    end
end

if(polyord>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(ind) = yin(i)*yin(j)*yin(k)*yin(l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyord>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(ind) = yin(i)*yin(j)*yin(k)*yin(l)*yin(m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

thresh = 2e-7;

if(laurentord>=1)
    for i=1:nVars
        if yin(i)>=thresh
            yout(ind) = 1/yin(i);
        else
            yout(ind) = 0;
        end;
        ind = ind+1
    end
end
if(laurentord>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            if yin(i)>=thresh && yin(j)>=thresh
                yout(ind) = 1/yin(i)/yin(j);
            else
                yout(ind) = 0;
            end;
            ind = ind+1;

        end
    end
end
if(laurentord>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                if yin(i)>=thresh && yin(j)>=thresh && yin(k)>=thresh
                    yout(ind) = 1/yin(i)/yin(j)/yin(k);
                else
                    yout(ind) = 0;
                end;
                ind = ind+1;
            end
        end
    end
end

if(absshifted>=1)
    for i=1:nVars
        yout(ind)=abs(yin(i)-1);
        ind = ind+1;
        yout(ind)=abs(yin(i)+1);
        ind = ind+1;
    end
end

if(sineord>=1)
    for k=1:sineord
        for i=1:nVars
            yout(ind) = sin(k*yin(i)) 
            ind = ind+1;
            yout(ind) = cos(k*yin(i));
            ind = ind+1;
        end
    end
end

if(intimed==1)
    yout(ind) = tin;
    ind = ind+1;
end

if(sintimed==1)
    yout(ind) = sin(tin);
    ind = ind+1;
    yout(ind) = cos(tin);
    ind = ind+1;
end


xprime=(yout*Xid)';


