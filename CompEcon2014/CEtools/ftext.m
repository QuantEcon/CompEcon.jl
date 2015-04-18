%% FTEXT
%
%  Inserts text in figure using a simplified format
%
%  Usage
%    ftext(x,y,txt,ha,va,fs)
%  Input
%    x     : x coordinate
%    y     : y coordinate
%    txt   : text to be inserted
%    ha    : horizontal allignment
%    va    : vertical allignment
%    fs    : fongsize of bullet (def: 18)

function ftext(x,y,txt,ha,va,fs)

if nargin<6, fs=18; end
text(x,y,txt,'HorizontalAlignment',ha,'VerticalAlignment',va,'FontSize',fs)