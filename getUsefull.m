function ret = getUsefull(iq, Su, Gi, f1stOfWhat)
    SpSym=Su+Gi;
    if strcmp(f1stOfWhat,'gi') == 1
        iq=iq(Gi+1:end);
    elseif strcmp(f1stOfWhat,'usefull')~=1 error('wrong f1stOfWhat'); end  
    i_end=(length(iq)+Gi)/SpSym;
    if floor(i_end)~=ceil(i_end) error("iq is not multiple of SpSym"); end
    
    res = [];
    for i=1:i_end
        res(i,:) = iq(1+(i-1)*SpSym:(i-1)*SpSym+Su);       
    end
    ret = res;
end