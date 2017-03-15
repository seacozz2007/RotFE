function s = ws2struct

vars = evalin('base', 'who')';
vals = evalin('base', ['{' sprintf('%s ', vars{:}) '}']);
c = [vars; vals];
s = struct(c{:});