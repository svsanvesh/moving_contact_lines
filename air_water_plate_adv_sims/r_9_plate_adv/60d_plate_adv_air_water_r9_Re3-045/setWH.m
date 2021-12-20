function  setWH(h, width,height )


pos = get(h, 'Position');
set(h, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size


% Here we preserve the size of the image when we save it.
set(h,'InvertHardcopy','on');
set(h,'PaperUnits', 'inches');
papersize = get(h, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(h,'PaperPosition', myfiguresize);


end

