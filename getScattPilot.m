function ret = getScattPilots(l)
    ret = [];
    for p=0:568
        k=0+3*mod(l,4)+12*p;
        if k > 6816 break; end
        ret = [ret k];
    end

end