function [map,lo,hi,prm] = cubehelix_view(N,start,rots,satn,gamma,irange,domain)
% An interactive figure for CubeHelix colormap parameter selection. With demo!
%
% (c) 2013-2020 Stephen Cobeldick
%
% View Dave Green's CubeHelix colorschemes in a figure.
%
% * Two colorbars give the colorscheme in color and grayscale.
% * A button toggles between 3D-cube and 2D-lineplot of the RGB values.
% * A button toggles an endless demo of randomly generated colorschemes.
% * Nine sliders allow real-time adjustment of the CubeHelix parameters.
% * Warning text if any RGB values are clipped.
% * Warning text if the grayscale is not monotonic increasing/decreasing.
%
%%% Syntax:
%  cubehelix_view
%  cubehelix_view(N)
%  cubehelix_view(N,start,rots,satn,gamma)
%  cubehelix_view(N,start,rots,satn,gamma,irange)
%  cubehelix_view(N,start,rots,satn,gamma,irange,domain)
%  cubehelix_view(N,[start,rots,satn,gamma],...)
%  cubehelix_view([],...)
%  cubehelix_view(axes/figure handles,...) % see "Adjust Colormaps"
%  [map,lo,hi] = cubehelix_view(...)
%
% Calling the function with an output argument blocks MATLAB execution until
% the figure is deleted: the final colormap and parameters are then returned.
%
% CubeHelix is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf
% For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
%
% Note: The original specification (the links above) misnamed the saturation
% option as "hue". In this function the saturation option is named "satn".
%
%% Adjust Colormaps of Figures or Axes %%
%
% Only R2014b or later. Provide axes or figure handles as the first input
% and their colormaps will be updated in real-time by CUBEHELIX_VIEW.
%
%%% Example:
%
% >> S = load('spine');
% >> image(S.X)
% >> cubehelix_view(gca)
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
%  N     = NumericScalar, an integer to define the colormap length.
%        = *[], colormap length of one hundred and twenty-eight (128).
%        = Array of axes/figure handles. R2014b or later only.
%  start = NumericScalar, the helix's start color, with R=1, G=2, B=3 (modulus 3).
%  rots  = NumericScalar, the number of R->G->B rotations over the scheme length.
%  satn  = NumericScalar, saturation controls how saturated the colors are.
%  gamma = NumericScalar, gamma can be used to emphasize low or high intensity values.
%  irange = NumericVector, range of brightness levels of the colormap's endnodes. Size 1x2.
%  domain = NumericVector, domain of the cubehelix calculation (endnode positions). Size 1x2.
%
%%% Outputs (these block code execution until the figure is closed!):
% map = NumericMatrix, the cubehelix colormap defined when the figure is closed.
% lo  = LogicalMatrix, true where <map> values<0 were clipped to 0. Size Nx3
% hi  = LogicalMatrix, true where <map> values>1 were clipped to 1. Size Nx3
%
% See also CUBEHELIX BREWERMAP RGBPLOT COLORMAP COLORMAPEDITOR COLORBAR UICONTROL ADDLISTENER

%% Input Wrangling %%
%
persistent ax2D ln2D ax3D pt3D txtH is2D cbAx cbIm pTxt pSld pVec
%
new = isempty(ax2D)||~ishghandle(ax2D);
dfn = 128;
upd = false;
upb = false;
hgv = [];
nmr = dfn;
%
err = 'First input must be a real positive scalar numeric or [].';
if nargin==0 || isnumeric(N)&&isequal(N,[])
	N = dfn;
elseif isnumeric(N)
	assert(isscalar(N),'SC:brewermap_view:NotScalarNumeric',err)
	assert(isnan(N)||isreal(N)&&isfinite(N)&&fix(N)==N&&N>0,...
		'SC:brewermap_view:NotRealPositiveNotNaN',err)
	N = double(N);
elseif all(ishghandle(N(:))) % R2014b or later
	assert(isgraphics(N(:),'axes')|isgraphics(N(:),'figure'),...
		'SC:cubehelix_view:NotAxesNorFigureHandles',...
		'First input may be an array of figure or axes handles.')
	hgv = N(:);
	nmr = arrayfun(@(h)size(colormap(h),1),hgv);
	N = nmr(1);
else
	error('SC:cubehelix_view:UnsupportedInput',err)
end
%
idp = 'SC:cubehelix_view:NotParameterVector';
idi = 'SC:cubehelix_view:NotLightRangeVector';
idd = 'SC:cubehelix_view:NotNodeDomainVector';
txp = '%s input can be a vector of the four CubeHelix parameters.';
txi = '%s input can be a vector of the endnode brightness levels (irange).';
txd = '%s input can be a vector of the endnode relative positions (domain).';
%
% Default pseudo-random parameters:
if nargin==0 || new
	clk = sum(clock*100);
	foo = @(n)min(3,max(0,sqrt(-log(rem(clk/pow2(n),1))*2)));
	bar = @(n)min(3,max(0,rem(clk/pow2(n),1)^36));
	pVec = [rem(clk,3);... % start
		min(3,max(-3,log10(rem(clk,1)/(1-rem(clk,1)))));... % rotations
		foo(5); foo(4);... % saturation and gamma
		bar(3); 1-bar(2); bar(1); 1-bar(0)]; % irange and domain
	% Original default parameters:
	% pVec = [0.5; -1.5;   1;   1; 0; 1; 0; 1];
	% %      [sta; rots; sat; gam; irng; domn]
end
% Parse input parameters:
switch nargin
	case 2
		pVec(1:4) = chvChk(4,start,idp,txp,'Second');
	case 3
		pVec(1:4) = chvChk(4,start,idp,txp,'Second');
		pVec(5:6) = chvChk(2,rots, idi,txi,'Third');
	case 4
		pVec(1:4) = chvChk(4,start,idp,txp,'Second');
		pVec(5:6) = chvChk(2,rots, idi,txi,'Third');
		pVec(7:8) = chvChk(2,satn, idd,txd,'Fourth');
	case 5
		pVec(1:4) = chvC2V(start,rots,satn,gamma);
	case 6
		pVec(1:4) = chvC2V(start,rots,satn,gamma);
		pVec(5:6) = chvChk(2,irange,idi,txi,'Sixth');
	case 7
		pVec(1:4) = chvC2V(start,rots,satn,gamma);
		pVec(5:6) = chvChk(2,irange,idi,txi,'Sixth');
		pVec(7:8) = chvChk(2,domain,idd,txd,'Seventh');
end
%
%% Ensure Figure Exists %%
%
% LHS and RHS slider bounds/limits, and slider step sizes:
lbd = [  1; 0;-3; 0; 0; 0; 0; 0; 0]; % left limit
rbd = [dfn; 3; 3; 3; 3; 1; 1; 1; 1]; % right limit
mnr = [100; 5; 5; 5; 5; 1; 1; 1; 1]; % minor step
mjr = [100; 5; 5; 5; 5; 1; 1; 1; 1]; % major step
%     [  N;st;ro;sa;ga;i1;i2;d1;d2]
stp = [mnr/100,mjr/10]; % [minor,major] step
%
% Define the 3D cube axis order:
xyz = 'RGB'; % choose order
[~,xyz] = ismember(xyz,'RGB');
%
% Parameter names for each slider:
spn = {'N';'start hue';'rotations';'saturation';'gamma';'irange(1)';'irange(2)';'domain(1)';'domain(2)'};
%
if new % Create a new figure.
	%
	% Figure parameters:
	M = numel(spn); % number of sliders
	gap = 0.01; % gaps
	bth = 0.04; % demo button height
	btw = 0.10; % demo button width
	uih = 0.40; % UI group height
	cbw = 0.23; % colorbar width (both together)
	cbh = 1-3*gap-bth; % colorbar height
	axh = 1-uih-2*gap; % axes height
	axw = 1-cbw-2*gap; % axes width
	slh = uih/M - gap; % slider height
	%
	figH = figure('HandleVisibility','callback', 'Color','white',...
		'IntegerHandle','off', 'NumberTitle','off',...
		'Name','CubeHelix Interactive Parameter Selector',...
		'MenuBar','figure', 'Toolbar','none', 'Tag',mfilename);
	%
	% Add 2D lineplot:
	ax2D = axes('Parent',figH, 'Position',[gap,uih+gap,axw,axh], 'Box','on',...
		'ColorOrder',[1,0,0; 0,1,0; 0,0,1; 0.6,0.6,0.6], 'HitTest','off',...
		'Visible','off', 'XLim',[0,1], 'YLim',[0,1], 'XTick',[], 'YTick',[]);
	ln2D = line([0,0,0,0;1,1,1,1],[0,0,0,0;1,1,1,1], 'Parent',ax2D,...
		'Visible','off', 'Linestyle','-', 'Marker','.');
	%
	% Add 3D scatterplot:
	ax3D = axes('Parent',figH, 'OuterPosition',[0,uih,axw+2*gap,1-uih],...
		'Visible','on', 'XLim',[0,1], 'YLim',[0,1], 'ZLim',[0,1], 'HitTest','on');
	pt3D = patch('Parent',ax3D, 'XData',[0;1], 'YData',[0;1], 'ZData',[0;1],...
		'Visible','on', 'LineStyle','none', 'FaceColor','none', 'MarkerEdgeColor','none',...
		'Marker','o', 'MarkerFaceColor','flat', 'MarkerSize',10, 'FaceVertexCData',[1,1,0;1,0,1]);
	view(ax3D,3);
	grid(ax3D,'on')
	lbl = {'Red','Green','Blue'};
	xlabel(ax3D,lbl{xyz(1)})
	ylabel(ax3D,lbl{xyz(2)})
	zlabel(ax3D,lbl{xyz(3)})
	%
	% Add warning text:
	txtH = text('Parent',ax2D, 'Units','normalized', 'Position',[0,1],...
		'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','r');
	%
	% Add demo button:
	demo = uicontrol(figH, 'Style','togglebutton', 'Units','normalized',...
		'Position',[1-cbw/2,1-bth-gap,cbw/2-gap,bth], 'String','Demo',...
		'Max',1, 'Min',0, 'Callback',@chvDemo); %#ok<NASGU>
	% Add 2D/3D button:
	is2D = uicontrol(figH, 'Style','togglebutton', 'Units','normalized',...
		'Position',[1-cbw/1,1-bth-gap,cbw/2-gap,bth], 'String','2D / 3D',...
		'Max',1, 'Min',0, 'Callback',@chv2D3D);
	%
	% Add colorbars:
	C(1,1,:) = [1,1,1];
	cbAx(2) = axes('Parent',figH, 'Visible','off', 'Units','normalized',...
		'Position',[1-cbw/2,gap,cbw/2-gap,cbh], 'YLim',[0.5,1.5],...
		'YDir','reverse', 'HitTest','off');
	cbAx(1) = axes('Parent',figH, 'Visible','off', 'Units','normalized',...
		'Position',[1-cbw/1,gap,cbw/2-gap,cbh], 'YLim',[0.5,1.5],...
		'YDir','reverse', 'HitTest','off');
	cbIm(2) = image('Parent',cbAx(2), 'CData',C);
	cbIm(1) = image('Parent',cbAx(1), 'CData',C);
	%
	% Add parameter sliders, listeners, and corresponding text:
	sv = mean([lbd,rbd],2);
	for m = M:-1:1
		Y = gap+(M-m)*(slh+gap);
		pLab(m) = uicontrol(figH,'Style','text', 'Units','normalized',...
			'Position',[gap,Y,btw,slh], 'String',spn{m}); %#ok<NASGU>
		pTxt(m) = uicontrol(figH,'Style','text', 'Units','normalized',...
			'Position',[gap+btw,Y,btw,slh], 'String','X');
		pSld(m) = uicontrol(figH,'Style','slider', 'Units','normalized',...
			'Position',[2*btw+gap,Y,axw-2*btw,slh], 'Min',lbd(m), 'Max',rbd(m),...
			'SliderStep',stp(m,:)/(rbd(m)-lbd(m)), 'Value',sv(m));
		addlistener(pSld(m), 'Value', 'PostSet',@(~,~)chvSldr(m));
	end
	%
end
%
%% Nested Functions %%
%
	function chvUpDt()
		% Update all graphics objects.
		%
		% Get CubeHelix colormap and grayscale equivalent:
		[map,lo,hi] = cubehelix(N, pVec(1:4),pVec(5:6),pVec(7:8));
		mag = map*[0.298936;0.587043;0.114021];
		%
		% Update colorbar values:
		set(cbAx, 'YLim',[0,abs(N)+(N==0)]+0.5);
		set(cbIm(1), 'CData',permute(map,[1,3,2]))
		set(cbIm(2), 'CData', repmat(mag,[1,1,3]))
		%
		% Update 2D line / 3D patch values:
		if  get(is2D, 'Value')
			set(ln2D, 'XData',linspace(0,1,abs(N)));
			set(ln2D,{'YData'},num2cell([map,mag],1).');
		else
			set(pt3D,...
				'XData',map(:,xyz(1)),...
				'YData',map(:,xyz(2)),...
				'ZData',map(:,xyz(3)),...
				'FaceVertexCData',map)
		end
		%
		% Update warning text:
		mad = diff(mag);
		wrn = {' Not Monotonic';' Clipped'};
		set(txtH,'String',wrn([any(mad<=0)&&any(0<=mad);any(lo(:))||any(hi(:))]));
		%
		% Update parameter value text:
		set(pTxt(1), 'String',sprintf('%.0f',N));
		set(pTxt(2:end), {'String'},sprintfc('%.2f',pVec));
		%
		% Update external axes/figure:
		nmr(1) = N;
		for k = 1:numel(hgv)
			colormap(hgv(k),cubehelix(nmr(k), pVec(1:4),pVec(5:6),pVec(7:8)));
		end
		%
		drawnow()
	end
%
	function chv2D3D(h,~)
		% Switch between 2D-line and 3D-cube.
		%
		if get(h,'Value') % 2D
			set(ax3D, 'HitTest','off', 'Visible','off')
			set(ax2D, 'HitTest','on',  'Visible','on')
			set(pt3D, 'Visible','off')
			set(ln2D, 'Visible','on')
		else % 3D
			set(ax2D, 'HitTest','off', 'Visible','off')
			set(ax3D, 'HitTest','on',  'Visible','on')
			set(ln2D, 'Visible','off')
			set(pt3D, 'Visible','on')
		end
		%
		chvUpDt();
	end
%
	function chvSldr(m)
		% Get new slider value.
		%
		if ~upd
			return
		end
		%
		val = get(pSld(m),'Value');
		if m==1
			N = round(val);
		else
			pVec(m-1) = val;
		end
		%
		chvUpDt()
	end
%
	function chvDemo(h,~)
		% Display random CubeHelix schemes.
		%
		% Parameter value step length:
		pvs = 0.03;
		% Functions to randomly specify new parameter values:
		randfn(7:8) = {@()rand(1,1).^42,@()1-rand(1,1).^42};
		randfn(3:4) = {@()sqrt(-log(rand(1,1))*2)};
		randfn(1:2) = {@()3*rand(1,1),@()randn(1,1)};
		randfn(5:6) = randfn(7:8);
		%
		% Define initial target values:
		tgt = pVec;
		%
		while ishghandle(h)&&get(h,'Value')
			%
			% new random targets:
			idr = abs(tgt-pVec)<1e-4;
			tgt(idr) = round(100*cellfun(@(fn)fn(),randfn(idr)))/100;
			% move onto close targets:
			idt = abs(tgt-pVec)<=pvs;
			pVec(idt) = tgt(idt);
			% move towards far targets:
			idm = ~(idr|idt);
			pVec(idm) = pVec(idm) + pvs*sign(tgt(idm)-pVec(idm));
			%
			upb = (upb || N<=1) && N<dfn;
			N = N - 1 + 2*upb;
			%
			% Update slider position:
			upd = false;
			try %#ok<TRYNC>
				set(pSld, {'Value'},num2cell(max(lbd,min(rbd,[N;pVec]))));
				chvUpDt();
			end
			upd = true;
			%
			% Faster/slower:
			pause(0.07);
		end
		%
	end
%
%% Initialize the Figure %%
%
set(pSld, {'Value'},num2cell(max(lbd,min(rbd,[N;pVec]))));
upd = true;
chvUpDt()
%
if nargout
	waitfor(ax2D);
	prm = pVec;
else
	clear map
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cubehelix_view
function V = chvC2V(varargin)
% Check that all of the input variables are real scalar numerics.
eid = 'SC:cubehelix_view:NotRealScalarNumeric';
msg = 'Input CubeHelix parameters must be %s.';
assert(all(cellfun(@isnumeric,varargin)),eid,msg,'numeric')
assert(all(cellfun(@isscalar,varargin)),eid,msg,'scalar')
assert(all(cellfun(@isfinite,varargin)),eid,msg,'finite')
assert(all(cellfun(@isreal,varargin)),eid,msg,'real')
V = cellfun(@double,varargin(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%chvC2V
function x = chvChk(n,x,eid,msg,ord)
% Check that the input variable <x> is real numeric vector with <n> elements.
assert(isnumeric(x)&&isreal(x)&&all(isfinite(x))&&isvector(x)&&numel(x)==n,eid,msg,ord)
x = double(x(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%chvChk
%
% Copyright (c) 2013 Stephen Cobeldick
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license