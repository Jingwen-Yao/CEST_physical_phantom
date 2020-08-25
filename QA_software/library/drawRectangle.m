function [] = drawRectangle(index,color)

% index 4 x 2
% 1 -------- 2
% |          |
% 3 -------- 4
line([index(1,1) index(2,1)],[index(1,2) index(2,2)],...
    'Color',color,'LineWidth',2,'LineStyle','--');
line([index(1,1) index(3,1)],[index(1,2) index(3,2)],...
    'Color',color,'LineWidth',2,'LineStyle','--');
line([index(2,1) index(4,1)],[index(2,2) index(4,2)],...
    'Color',color,'LineWidth',2,'LineStyle','--');
line([index(3,1) index(4,1)],[index(3,2) index(4,2)],...
    'Color',color,'LineWidth',2,'LineStyle','--');

end