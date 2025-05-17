function map = makecmap(c,n,cspace,t0)
import c19n_subfunctions.*

% cspace can be 'rgb', 'lab' or 'xyz'.
if nargin < 3
	cspace = 'rgb';
end
m = size(c,1);
if nargin < 4
	t0 = linspace(0,1,m)';
else
	t0 = t0(:)/t0(end);
end
t = linspace(0,1,n)';
switch cspace
	case 'rgb'
		map = interp1(t0,c,t);
	case {'lab','xyz'}
		xform1 = makecform(['srgb2',cspace]);
		xform2 = makecform([cspace,'2srgb']);
		lab = applycform(c,xform1);
		map_other = interp1(t0,lab,t);
		map = applycform(map_other,xform2);
	otherwise
		error('Unknown color space, %s',cspace)
end
map = min(max(map,0),1);
end
