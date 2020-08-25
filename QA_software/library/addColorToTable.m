function [] = addColorToTable(z_score, column_number, uit)

z_score(z_score > 5) = 5;
z_score(z_score < -5) = -5;

for row = 1:length(z_score)
    if z_score(row) >= 0
        cellColor = [1 1 1] - [0 0.3 0.3]/5*z_score(row);
    else
        cellColor = [1 1 1] + [0.3 0.3 0]/5*z_score(row);
    end
    s = uistyle('BackgroundColor',cellColor);
    addStyle(uit,s,'cell',[row,column_number]);    
end

kU = find(z_score > 1.96);
kL = find(z_score < -1.96);
if ~isempty(kU)
    s = uistyle('FontColor','red');
    addStyle(uit,s,'cell',[kU,ones(size(kU))*column_number]);
end
if ~isempty(kL)
    s = uistyle('FontColor','blue');
    addStyle(uit,s,'cell',[kL,ones(size(kL))*column_number]);
end

end