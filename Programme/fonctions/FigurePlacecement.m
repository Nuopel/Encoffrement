function [ ] = FigurePlacecement(opt)
% [ ] = figurePlacecement( )
% Create the number of  figure reparted in the display 

n =  length(findobj('type','figure'));
%% window figure management
screen_size = get(0,'ScreenSize') ;
pcw=screen_size(3);
pch=screen_size(4);

%%
if nargin < 1
row  = ceil(sqrt(n));
line = floor(sqrt(n));
else
row  = floor(sqrt(n));
line = ceil(sqrt(n));
end

if row*line < n
    row=row+1;
end
%%
widthw=pcw/row;
heightw=pch/line;
for ii=1:line
    for jj=1:row
    f1=figure((ii-1)*row+jj);
    set(f1, 'Position', [widthw*(jj-1),heightw*(ii-1), widthw-7,heightw-88]);
    end
end
nf2 =  findobj('type','figure');
for ii=1:length(nf2)-n
    close(num2str(n+ii))
end
end

