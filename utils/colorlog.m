function colormap = colorlog(x)
%computes a log color, allowing for 0s and negative values

fact=1000;

if(x<0)
    colormap=-log10(-x*1000);
elseif(x==0)
    colormap=0;
else
    colormap=log10(x*1000);
end;

end
