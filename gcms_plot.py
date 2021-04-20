import pyteomics.mzxml
import pandas as pd
import numpy as np
import bokeh.io
import bokeh
import bokeh.plotting

from random import random
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, CustomJS, LinearAxis
from bokeh.models.tools import*
from bokeh.models.ranges import Range1d
from bokeh.models import CustomJS

from bokeh.plotting import figure, output_file, show

def gcms_plot(df_chromat, df_ms, plot_height=400, plot_width=800):

    x = df_chromat['num']
    y = df_chromat['chromat']
    t = df_chromat['t']

    s1 = ColumnDataSource(data=dict(x=x, y=y, t=t))


    TOOLS = "pan,box_zoom,wheel_zoom,save,reset,box_select"
    min_ind = x.astype(int).min()
    p1 = figure(plot_width=plot_width,
                plot_height=plot_height,
                x_range=(x.astype(int).min(), x.astype(int).max()),
                x_axis_label='Retention time',
                y_axis_label='Counts',
                title='Chromatogram (double click on plot for mass spec @ retention time)',
                tools=TOOLS)

    p1.line('x', 'y', source=s1, alpha=0.6)
    p1.circle('x', 'y', source=s1, alpha=0.6, size=1, line_alpha=0)
    p1.xaxis.visible = False
    p1.extra_x_ranges = {"Retention time": Range1d(start=t.min(), end=t.max())}
    p1.add_layout(LinearAxis(x_range_name="Retention time"), 'below')

    s4 = ColumnDataSource(data=dict(x=[], y=[]))
    p1.line('x', 'y', source=s4, color='black')

    s2 = ColumnDataSource(data=dict(x=[], y=[]))
    TOOLTIPS = [
        ("(x,y)", "($x, $y)"),
    ]
    p2 = figure(plot_width=plot_width, plot_height=plot_height, x_range=(0, 600), y_range=(0, 5000), title="Mass spec", tooltips=TOOLTIPS, x_axis_label='Mass', y_axis_label='Counts')
    s5 = ColumnDataSource(data=dict(x=[], y=[]))
    p2.vbar('x', 0.5, 'y', source=s2, alpha=0.6)
    p2.circle('x', 'y', source=s5, alpha=0.6)
    p2.js_on_event('reset', CustomJS(args=dict(s2=s2, yr=p2.y_range), code="""
            var d2 = s2.data;
            yr.start=Math.min(...d2['y']);
            yr.end=Math.max(...d2['y'])*1.01;
            yr.change.emit();
        """)
    )
    s3 = ColumnDataSource(data=dict(x=df_ms['mass'], y=df_ms['intensity'], i=df_ms['num']))


    s1.selected.js_on_change(
        "indices",
        CustomJS(
            args=dict(s1=s1, s2=s2, s3=s3, s5=s5, title=p2.title, min_ind=min_ind),
            code="""
            var inds = cb_obj.indices;
            var start = Math.min(...inds);
            var stop = Math.max(...inds);
            var d5 = s5.data;
            var d3 = s3.data;
            if (start > 0 && start < 10000000) {
                var x_new = [];
                var y_new = [];
                for (var i=start; i<stop+1; i++) {
                    x_new = x_new.concat(d3['x'][i]);
                    y_new = y_new.concat(d3['y'][i]);
                }
                var sum = 0;
                for (var i=0; i<y_new.length; i++) {
                    sum += y_new[i];
                }
                d5['x'] = x_new;
                d5['y'] = y_new;
                //title.text=" TIC: ".concat(sum.toString()).concat(" Peak locs: ").concat(start.toString()).concat(", ").concat(stop.toString());
            }
            s5.change.emit();
            title.change.emit()
        """,
        ),
    )

                #d2['x'] = [].concat(d3['x'].slice(start,stop+1))
                #d2['y'] = [].concat(d3['x'].slice(start,stop+1))

    p1.js_on_event(bokeh.events.DoubleTap, CustomJS(args=dict(s1=s1, s2=s2, s3=s3, s4=s4, s5=s5, yr=p2.y_range, title=p2.title, min_ind=min_ind), code="""
            var inds = cb_obj.indices;
            var d1 = s1.data;
            var d2 = s2.data;
            var d3 = s3.data;
            var d4 = s4.data;
            var d5 = s5.data;
            var i = Math.round(cb_obj.x);

            d4['x'] = [i, i];
            var y_max = Math.max(...d1['y']);
            d4['y'] = [0, y_max*1.01];
            
            i = i-min_ind;
            
            d2['x'] = d3['x'][i];
            d2['y'] = d3['y'][i];

            d5['x'] = d3['x'][i];
            d5['y'] = d3['y'][i];

            yr.start=Math.min(...d3['y'][i]);
            yr.end=Math.max(...d3['y'][i])*1.01;
            title.text="Mass spec, retention time: ".concat(d1['t'][i].toString());

            yr.change.emit();
            s2.change.emit();
            s4.change.emit();
            s5.change.emit();
            title.change.emit()
        """)
    )
    #Math.round(cb_obj.x)
    layout = column(p1, p2)
    return(layout)