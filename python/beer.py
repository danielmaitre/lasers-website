import numpy as np
from bokeh.layouts import layout, widgetbox
from bokeh.models import CustomJS, Slider
from bokeh.models.glyphs import Circle, Line
from bokeh.plotting import figure, output_file, show, ColumnDataSource

I_0 = 1.0 # intensity amplitude

''' Population density of the two levels in the laser transition in m^-3'''
N_1 = 1e16
N_2 = 10e16

''' Degeneracies of the two levels in the laser transition '''
g_1 = 1
g_2 = 1

NStar = N_2 - (g_2/g_1)*N_1 # population inversion density in m^-3
sigma = 4e-17 # in m^2, see specHoles_stripped.py and/or lecture notes, should get equation for it depending on omega

gran = 100
zs = np.linspace(0, 1.0, gran) # distance through gain medium, in m
Is = np.zeros((gran))
for i in range(gran):
    Is[i] = I_0 * np.exp(sigma * NStar * zs[i])



popScale = 1e16
''' callback data sources '''
attributes = ColumnDataSource(data=dict(sigma=[sigma], scale=[popScale]))
plotArrs = ColumnDataSource(data=dict(zs=zs, Is=Is))
flatArrs = ColumnDataSource(data=dict(x=[zs[0], zs[gran-1]], y=[I_0]*2))


''' ************** create plot **************** '''
maxI = np.max(Is)
                              
p1 = figure(title='Intensities in gain medium',
            plot_width=400, plot_height=400,
            y_range = (0, 3) # note widths different to holeBurning (change to 400,400 if display issues)
            )
p1.xaxis.axis_label = "Propagation distance (m)"
p1.yaxis.axis_label = "Transmitted Intensity in units of initial intensity, (I_T)/(I_0)"
p1.line('zs', 'Is', source=plotArrs, line_width=3.0, line_alpha=1.0)
p1.line('x', 'y', source=flatArrs, line_width=1.0, line_alpha=0.8, line_dash = 'dashed', line_color='black')
#p1.legend.location = "top_left"
#p1.legend.background_fill_alpha = 0.5

''' ****************************** '''

callback = CustomJS(args=dict(srce1=attributes, srce2=plotArrs, srce3=flatArrs), code = """
    var data1 = srce1.data;
    sigma = data1['sigma'][0];
    var popScale = data1['scale'][0];

    var data2 = srce2.data;
    zs = data2['zs'];
    Is = data2['Is'];
    gran = zs.length;

    var data3 = srce3.data;
    flatI = data3['y'];

    var N_1 = N1.value*popScale;
    var N_2 = N2.value*popScale;
    I_0 = 1.0;//var I_0 = I0.value;

    flatI[0] = I_0; flatI[1] = I_0;

    var NStar = N_2 - (g2.value/g1.value)*N_1;
    for (var i = 0; i < gran; i++) {
        Is[i] = I_0 * Math.exp(sigma * NStar * zs[i]);
    }
    
    srce2.change.emit();
    srce3.change.emit();
    
""")






maxPop = 10e16
maxDegeneracy = 10
N1Slider = Slider(start=0, end=maxPop/popScale, value=N_1/popScale, step=maxPop/(1000*popScale),
                      title='$N_1$ (population density of level 1, in units of 1e16 m$^{-3}$)',
                      callback=callback,
                      )
N2Slider = Slider(start=0, end=maxPop/popScale, value=N_2/popScale, step=maxPop/(1000*popScale),
                    title='$N_2$ (population density of level 2, in units of 1e16 m$^{-3}$)',
                    callback=callback,
                    )
g1Slider = Slider(start=1, end=maxDegeneracy, value=g_1, step=1,
                    title='$g_1$ (degeneracy of level 1)',
                    callback=callback,
                    )
g2Slider = Slider(start=1, end=maxDegeneracy, value=g_2, step=1,
                   title = '$g_2$ (degeneracy of level 2)',
                   callback=callback)
#I0Slider = Slider(start=1.0, end=15.0, value=I_0, step=0.01, # add a modulus so a multiple of step?
                   #title = '$I_0$, intensity in absence of laser field, in W m$^{-2}$',
                   #callback=callback)


callback.args['N1'] = N1Slider
callback.args['N2'] = N2Slider
callback.args['g1'] = g1Slider
callback.args['g2'] = g2Slider
#callback.args['I0'] = I0Slider

'''
div = Div(text="""<div id="FSR_info">
    Free spectral range $\Delta\omega_{{FSR}}$: {0:3g} GHz
</div>""".format(deltaOMEGA_FSR*unitScale), width=400, height = 20)
'''

l = layout([[p1],
  #[p1, p2],
  [widgetbox(N1Slider, N2Slider, g1Slider, g2Slider
             #, I0Slider
             )],
], sizing_mode='stretch_both')
output_file("beer.html", title="Beer Law")
show(l)
