function Num = field2num(Field,DefaultNum)
% Num = field2num(Field [,DefaultNum]) reads numerical input in Field (scalar or
% a two-element vector separated by colon) into Num relayed in the field. If input
% is invalid and DefaultNum specified, Num and the field are set to DefaultNum.

% if ~isequal(Field,'hObject'); eval('Field=[''handles.'' Field];'); end % does not work
Str=get(Field,'String'); if ~isempty(findstr(Str,' ')); Str='XXX'; end; Str=strrep(Str,':',' '); 
% Num=sort(str2num(Str));
Num=str2num(Str);
if nargin>1
   if isempty(Num)||length(Num)>2; Num=DefaultNum;
      if length(DefaultNum)==1; set(Field,'String',sprintf('%g',DefaultNum))
         else set(Field,'String',sprintf('%g:%g',DefaultNum))
      end
   end
end