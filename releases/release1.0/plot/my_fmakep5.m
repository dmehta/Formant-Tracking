% ======================================================================
%                            FMAKEPLOT.M
% ======================================================================
%
%      This routine makes an existing plot into a viewgraph.
%
%   PROGRAMMER  : S. L. Wilson
%
%   DATE   CODE : 02 December  1993
%   UPDATE CODE : 15 December  1994 - Modified to do multiple subplots
%                    and you can hit return at the VG subtitle
%                    input line and it moves the title down and
%                    doesn't write a subtitle.  If there are two
%                    axes, for instance a plot and a color scale,
%                    it squishes them vertically to fit the title
%                    and the logo.  If there are more than two
%                    axes, like four subplots, it squishes in both
%                    directions.  These are hard-wired fixes so
%                    they may have to be played with to get a
%                    truly pleasing result.
%                 17 May 1995 - Turned into a function so user can
%                    give the title, subtitle, and serial in a function
%                    call.
%
% ----------------------------------------------------------------------
%
%   USAGE  :  fmakeplot(vg_title,vg_subtitle)
%
% ======================================================================

%function [ a ]            =  fmakevg(vg_title,vg_subtitle)  ;


%
%---> Define text font
%
text_font = 'Arial';

%
%---> Define text sizes.
%
system_dependent(14,'on');
text_size                 =  12                                        ;
label_size                =  16                                        ;
vg_title_size             =  24                                        ;
vg_subtitle_size          =  20                                        ;
vg_serial_size            =  6                                         ;
mark_size                 =  12                                        ;

%
%---> Get the current figure handle and find number and handle to
%---> children of the figure i.e., the axes present in the figure.
%

h_fig                     =  gcf                                       ;
axes_vec                  =  get(h_fig,'Children')                     ;
%figure_window             =  [ 0.0 0.0 1.0 1.0 ]                       ;
%htitle_axes               =  axes('Position',figure_window,          ...
%                                  'Visible','Off')                     ;

%
%---> Loop over figure's children (the axes).
%

for kk = 1 : length(axes_vec)

    h_figaxes             =  axes_vec(kk)  ;                            

    if ( length(axes_vec) == 1 )

%       figure_window      =  [ 0.20 0.20 0.60 0.50 ]                   ;
%       figure_window      =  [ 0.30 0.20 0.40 0.50 ] ; % This makes figure square.

    elseif ( length(axes_vec) == 2 )
%
%---> If bottom is over 1/2 lower it. If under 1/2 raise it.
%
       rect_loc           =  get(h_figaxes,'Position')                 ;

%       figure_window      =  rect_loc                                  ;

    else

%
%---> If bottom is over 1/2 lower it. If under 1/2 raise it a lot.
%
       rect_loc           =  get(h_figaxes,'Position')                 ;
%       figure_window      =  rect_loc                                  ;

    end

%
%---> Get the handles to the children of the current axes.
%---> These will be of type 'text' and 'line'.
%---> Make the lines thick, 2 points, and Markers big.
%---> Make the text bold and big.
%
    h_children            =  get(h_figaxes,'Children')                 ;
    n_children            =  length(h_children)                        ;

	
    for ii = 1 : n_children
        child_type        =  get(h_children(ii),'Type')                ;
        if     ( strcmp(child_type,'Line') == 1 |                    ...
                 strcmp(child_type,'line') == 1 |                    ...
                 strcmp(child_type,'LINE') == 1 )
           set(h_children(ii),'LineWidth',2,                         ...
                              'ButtonDownFcn','setlinec',            ...
                              'MarkerSize',9)                          ;
        elseif ( strcmp(child_type,'Text') == 1 |                    ...
                 strcmp(child_type,'text') == 1 |                    ...
                 strcmp(child_type,'TEXT') == 1 )
	   if (strcmp(get(h_children(ii),'FontUnits'),'normalized'))
		set(h_children(ii),'FontUnits','points');
	   end
           set(h_children(ii),'FontName',text_font,                ...
                              'FontWeight','Normal',                   ...
                              'FontSize',text_size)                    ;
        end
    end
	if  strcmp(get(h_figaxes,'Tag'),'legend') == 0
%
%---> Make the axes numbers bold and set the position to figure_window.
%---> Box should be 'On' or 'Off' as the user prefers.
%

if ~strcmp(get(h_figaxes,'Type'),'uicontrol')
    set(h_figaxes,'FontName',text_font,                            ...
                  'FontSize',text_size,                              ...
                  'FontWeight','Normal',                               ...
                  'LineWidth',2,                                     ...
                   'Box','On')                                          ;
%                  'Position',figure_window,                          ...

%
%---> Make all labels big and bold.
%
    set(get(h_figaxes,'Xlabel'),'FontName',text_font,              ...
                                'FontSize',label_size,               ...
                                'FontWeight','Normal')                   ;

    set(get(h_figaxes,'Ylabel'),'FontName',text_font,              ...
                                'FontSize',label_size,               ...
                                'FontWeight','Normal')                   ;

    set(get(h_figaxes,'Zlabel'),'FontName',text_font,              ...
                                'FontSize',label_size,               ...
                                'FontWeight','Normal')                   ;

    set(get(h_figaxes,'Title'), 'FontName',text_font,              ...
                                'FontSize',label_size,               ...
                                'FontWeight','Normal')                   ;

end
end
end

%end

%
%---> Make sure to set orientation to landscape.
%

%orient landscape

%
%---> C'est fini!
%
