function [yout] = poolData3(tin,yin,nVars,inoned,polyorder,...
    laurentord,absshifted,sineord,intimed,sintimed)

% Copyright 2015 to 2019, All Rights Reserved
% Code initially by Steven L. Brunton
% for Paper, "Discovering Governing Equations from Data:
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% Modifications Robert Niven from June 2018 to Jan 2024, to allow for
% - better handle polynomials, laurent functions and sine functions
% - include other alphabet choices including shifted absolute values
% - change plot label output


m = size(yin,1);

ind = 1;

% poly order 0
if(polyorder>=0)
    if (inoned==1)
         yout(:,ind) = ones(m,1);
         %funstr{ind,1}  = '1';
         ind = ind+1;
    end
end
% poly order 1
if(polyorder>=1)
    for i=1:nVars
        yout(:,ind) = yin(:,i);
        ind = ind+1;
    end
end
% poly order 2
if(polyorder>=2)
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;         
        end
    end    
end
% poly order 3
if(polyorder>=3)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end  
end
% poly order 4
if(polyorder>=4)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end
% poly order 5
if(polyorder>=5)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
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
        ithresh = yin(:,i)>=thresh;
        yout(ithresh,ind) = 1./yin(ithresh,i);
        yout(~ithresh,ind) = 0;
        %yout(yin(:,i)>=thresh,ind) = 1./yin(yin(:,i)>=thresh, i);
        %yout(yin(:,i)<thresh, ind) = 0;
        ind = ind+1
    end
end
% laur order 2
if(laurentord>=2)
    for i=1:nVars
        for j=i:nVars
            ijthresh = yin(:,i)>=thresh & yin(:,j)>=thresh;
            yout(ijthresh,ind) = 1./yin(ijthresh,i)./yin(ijthresh,j);
            yout(~ijthresh,ind) = 0;
            ind = ind+1;
            
        end
    end
end
% laur order 3
if(laurentord>=3)
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                ijkthresh = yin(:,i)>=thresh & yin(:,j)>=thresh & yin(:,k)>=thresh;
                yout(ijkthresh,ind) = 1./yin(ijkthresh,i)./yin(ijkthresh,j)./yin(ijkthresh,k);
                yout(~ijkthresh, ind)= 0;
                ind = ind+1; 
            end
        end
    end           
end
    
if(laurentord==1) && (polyorder==2)
    %special case for Watt system; have been careful with combinatorics
    for i=1:nVars
        for j=i:nVars
            for k=1:nVars
                 if (k~=i) && (k~=j)
                    kthresh = yin(:,k)>=thresh;
                    yout(kthresh,ind) = yin(kthresh,i).*yin(kthresh,j)./yin(kthresh,k);
                    yout(~kthresh,ind)= 0;
                    ind = ind+1;
                 end
            end
        end
    end
end

if(absshifted>=1)
    for i=1:nVars
        yout = [yout abs(yin(:,i)-1) abs(yin(:,i)+1)];
    end
end

if(sineord>=1)
    for k=1:sineord
        for i=1:nVars
            yout = [yout sin(k*yin(:,i)) cos(k*yin(:,i))];
        end
    end
end

if(intimed==1)
    yout = [yout tin];
end

if(sintimed==1)
    yout = [yout sin(tin) cos(tin)];
end

