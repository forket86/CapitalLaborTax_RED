function  temp_output=SimNetwork(X,ymaxin,yminin,xmaxin,xminin,ymaxout,yminout,xmaxout,xminout,IW,B1,LW,b2)
    in = (ymaxin-yminin)*(X-xminin)./(xmaxin-xminin) + yminin;
    h = tansig(IW*in+B1);
    y2n = LW*h + b2;
    temp_output=(xmaxout-xminout).*(y2n-yminout)./(ymaxout-yminout) + xminout;
end