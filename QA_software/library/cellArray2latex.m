function [latex] = cellArray2latex(cellArray, input, align, fixWidth)

if nargin < 3
    align = 'c';
end

if nargin < 4
    header = ['\begin{tabular}','{',repmat([align '|'],1,size(cellArray,2)-1),align,'}'];
else
    header = ['\begin{tabular}',fixWidth];
end

latex = {header};
hLine = ['\hhline{',repmat('~|',1,size(cellArray,2)-1),'~}'];

% generate table
if input.booktabs
    latex(end+1) = {'\toprule'};
end
for i = 1:size(cellArray,1)
    if i == 2 && input.booktabs
        latex(end+1) = {'\midrule'};
    end
    if input.tableBorders
        latex(end+1) = {hLine};
    end
    rowStr = '';
    for j = 1:size(cellArray,2)
        dataValue = cellArray{i,j};
        if iscell(dataValue)
            dataValue = dataValue{:};
        elseif isnan(dataValue)
            dataValue = input.dataNanString;
        elseif isnumeric(dataValue)
            dataValue = num2str(dataValue,dataFormatArray{i,j});
        end
        if j == 1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    latex(end+1) = {[rowStr,' \\']};
end

if input.booktabs
    latex(end+1) = {'\bottomrule'};
end

if input.tableBorders
    latex(end+1) = {hLine};
end
latex(end+1) = {'\end{tabular}'};

latex = cellfun(@(x) strrep(x,'#','\#'),latex,'UniformOutput',false);
latex = cellfun(@(x) strrep(x,'±','$\pm$'),latex,'UniformOutput',false);
latex = cellfun(@(x) strrep(x,'%','\%'),latex,'UniformOutput',false);

end

