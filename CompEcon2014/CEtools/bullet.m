%% BULLET
%
%  Inserts centered bullet in figure using simplified format
%
%  Usage
%    bullet(x,y,size,color)
%  Input
%    x     : x coordinate
%    y     : y coordinate
%    size  : size of bullet (def: 25)
%    color : color of bullet (def: black)

function bullet(x,y,size,color)

if nargin<3, size=25; end
if nargin<4, color='k'; end
text(x,y,'\bullet','color',color,'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',size)