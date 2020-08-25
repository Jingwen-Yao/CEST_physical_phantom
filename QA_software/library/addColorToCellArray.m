function [Tplot] = addColorToCellArray(z_score, column_number, Tplot, zcutoff)

z_score(z_score > 4) = 4;
z_score(z_score < -4) = -4;

for row = 1:length(z_score)
    if z_score(row) >= 0
        cellColor = ['\cellcolor{red!' num2str(ceil(25/4*z_score(row)),'%d') '}'];
    else
        cellColor = ['\cellcolor{blue!' num2str(ceil(-25/4*z_score(row)),'%d') '}'];
    end
    if z_score(row) > zcutoff
        textColor = '\textcolor{red}';
    elseif z_score(row) < -zcutoff
        textColor = '\textcolor{blue}';
    else
        textColor = '';
    end
    Tplot{row+1,column_number} = ...
        [textColor '{' cellColor '{' Tplot{row+1,column_number} '}}'];
end

end