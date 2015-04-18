function printfigures(nam,fignum,BoxOff,TitleOff)
% Save figures as .eps files for overheads and manuscripts
% Copyright (c) 1997-2014, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

% for i=1:fignum
%   figure(i)
%   legend('boxoff')
% end
% return

% return

if nargin<3, BoxOff = 1; end
if nargin<4, TitleOff = 1; end
for i=1:fignum
  figure(i)
  legend('boxoff')
  if BoxOff, box off; end
  if TitleOff, title([]); end
  eval(['print -depsc ' nam num2str(i)])
end