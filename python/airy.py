import numpy as np
import math
from bokeh.layouts import layout, row, column, widgetbox
from bokeh.models import CustomJS, Slider, Dropdown, Div
from bokeh.models.glyphs import Circle, Line, CircleCross
from bokeh.plotting import figure, output_file, show, ColumnDataSource#, TapTool

pi = np.pi
c = 3e8
M_A = 40 # relative atomic mass of Ar
lambd = 488e-9 # laser wavelength

L = 0.2 # length of cavity in m
T = 1200 + 273 # discharge T in kelvin (check whether function of atom or not)
R_1 = 0.95 # reflectivity of mirrors
R_2 = 0.95
offset = 0; # frequency offset of zeroth mode from zero detuning
eta = 1.0 # refractive index of gas in etalon

''' plotting variables '''
gran = 100000 # how detailed to make plot (reduce to increase speed)
unitScale = 1/(2*pi*1e9) # Ghz
plotWidth = 1.5 # how much of gaussian to show

R = (R_1*R_2)**0.5
omega_12 = 2*pi*c / lambd
deltaOMEGA_I = 7.16e-7 * omega_12 * (T / M_A)**(0.5)
detunings = np.linspace(-deltaOMEGA_I*plotWidth, deltaOMEGA_I*plotWidth, gran)
detGHz = detunings * unitScale # scale to plotting units

deltaOMEGA_FSR = pi*c/(eta*L)
finesse = pi * np.sqrt(R) / (1-R)
moddedOffset = offset % deltaOMEGA_FSR # calculate closest mode peak to right of resonance

def airy(det):
    return 1/(1+(4*finesse**2/pi**2)*(np.sin(det/2))**2)


deltas = np.linspace(-2*pi*deltaOMEGA_I*plotWidth/deltaOMEGA_FSR, 2*pi*deltaOMEGA_I*plotWidth/deltaOMEGA_FSR, gran)
deltas -= 2*pi*moddedOffset/deltaOMEGA_FSR
As = airy(deltas) # Airy function due to cavity width
maxA = np.max(As)
As /= maxA # scale and position on plot

refStep = 0.01
minRef = refStep
maxRef = 0.99
possibleRefs = np.arange(minRef, maxRef+refStep, refStep)
Rs = np.array(())
Fs = np.array(())

RLoc = 0
i=0
#print(R)
for R_i in possibleRefs:
    for R_j in possibleRefs:
        #print('{}, {}'.format(R_i, R_j))
        currR = (R_i*R_j)**0.5
        Rs = np.append(Rs, currR)
        Fs = np.append(Fs, pi*currR**0.5 / (1 - currR)) # finesse
        if (np.isclose(R_1, R_i)) & (np.isclose(R_2, R_j)):
            RLoc = i
            #print(RLoc)
        i+=1
            
#print(RLoc)

markerSizeUnselect = 3
markerSizeSelect = 15
sizes = [markerSizeUnselect]*Fs.size
sizes[RLoc] = markerSizeSelect

''' callback data sources '''
attributes = ColumnDataSource(data=dict(off=[offset], sSize=[markerSizeSelect], uSize=[markerSizeUnselect],
                                        uScale=[unitScale]))
plot1Arrs = ColumnDataSource(data=dict(Rs=Rs, Fs=Fs, sizes=sizes))
plot2Arrs = ColumnDataSource(data=dict(det=detunings, As=As, detGHz=detGHz))

''' ************** create plot **************** '''
lineSelection = Line(line_color='black', line_width=1, line_alpha=1.0)
lineUnselection = Line(line_color='black', line_width=1, line_alpha=1.0)
circleSelection = CircleCross(fill_color='red', fill_alpha=0.3, line_color='red', line_alpha=0.5, line_width=1,
                              angle=45)
circleUnselection = CircleCross(fill_alpha=0.1, fill_color='blue', line_color=None)

p1 = figure(title='Finesse of etalon against reflectivity (red cross indicates current finesse value)',
            name = 'finesseScatter',
            plot_width=300, plot_height=400) # note widths different to holeBurning (change to 400,400 if display issues)
p1.xaxis.axis_label = "Overall reflectivity of etalon"
p1.yaxis.axis_label = "Coefficient of finesse"
r1=p1.circle_cross('Rs', 'Fs', source=plot1Arrs, line_color=None, fill_color='blue', size='sizes', fill_alpha=0.1)
r1.selection_glyph = circleSelection
r1.nonselection_glyph = circleUnselection

p2 = figure(x_range=(-deltaOMEGA_I*plotWidth*unitScale,deltaOMEGA_I*plotWidth*unitScale),
            #y_range=(0,1),
            title='Normalised Airy function of a Fabry-Perot etalon',
            plot_width=500, plot_height=400) # note widths different to holeBurning (change to 400,400 if display issues)
p2.xaxis.axis_label = "detuning from laser transition frequency, in GHz"
p2.yaxis.axis_label = "Normalised transmitted intensity"
r2=p2.line('detGHz', 'As',
        source=plot2Arrs,
        line_width=1, line_alpha=1.0, color='black', name = 'airyLine')
r2.selection_glyph = lineSelection
r2.nonselection_glyph = lineUnselection
''' ****************************** '''

callback = CustomJS(args=dict(srce1=attributes, srce2=plot2Arrs, srce3=plot1Arrs), code = """
    var c = 3e8;

    var data1 = srce1.data;
    offset = data1['off'][0];
    var unselectSize = data1['uSize'][0];
    var selectSize = data1['sSize'][0];
    var unitScale = data1['uScale'][0];

    var data2 = srce2.data;
    detunings = data2['det'];
    As = data2['As'];
    gran = detunings.length;

    var data3 = srce3.data;
    Rs = data3['Rs'];
    Fs = data3['Fs'];
    sizes = data3['sizes'];
    
    var L = length.value;
    var eta = etaSl.value;
    
    var R_1 = ref1.value;
    var R_2 = ref2.value;
    var RStart = ref1.start;
    var REnd = ref1.end;
    var RStep = ref1.step;
    

    var numRefValues = Math.round((REnd - RStart) / RStep)+1;
    var R_1Idx = Math.round((R_1 - RStart) / RStep) * numRefValues;
    var R_2Idx = Math.round((R_2 - RStart) / RStep);
    var selectedIdx = R_1Idx + R_2Idx;
    var selectIndices = [selectedIdx];
    console.log(selectIndices);
    
    var mySelected = srce3.selected;
    mySelected['1d']['indices'] = selectIndices;
    srce3.selected = mySelected;

    for (var i=0; i<sizes.length; i++) {
        sizes[i] = unselectSize;
    }
    sizes[selectedIdx] = selectSize;
    
    srce3.change.emit();

    /* recalculate values */
    var deltaOMEGA_FSR = math.pi*c/(L*eta);
    var moddedOffset = offset % deltaOMEGA_FSR;
    R = Math.sqrt(R_1*R_2);
    var finesse = math.pi * Math.sqrt(R) / (1-R);

    /* recalculate the Airy function for the plot */
    newAs = new Array(gran);
    for (var i=0; i<gran; i++) {
        newAs[i] = airy(((2*math.pi/deltaOMEGA_FSR)*detunings[i] - 2*math.pi*moddedOffset/deltaOMEGA_FSR), finesse); 
    }
    maxA = math.max(newAs);
    for (var i=0; i<gran; i++) {
        As[i] = newAs[i]/maxA; // scale Airy function
    }
    
    srce2.change.emit();

    var fsrDiv = document.getElementById('FSR_info');
    fsrDiv.textContent = math.format(deltaOMEGA_FSR*unitScale, 6)
""")



'''

'''

airy = """
    function airy(det, finesse) {
        return 1/(1+(4*math.square(finesse)/(math.square(math.pi)))*math.square(math.sin(det/2)));
    }
"""



# two reflectivity sliders needed now
lengthSlider = Slider(start=0.01, end=0.3, value=L, step=0.01,
                      title='$L$ (cavity length in m)',
                      callback=callback,
                      #callback_policy='mouseup'
                      )
refl1Slider = Slider(start=minRef, end=maxRef, value=R, step=refStep,
                    title='$R_1$ (reflectivity of mirror 1)',
                    callback=callback,
                    #callback_policy = 'mouseup'
                    )
refl2Slider = Slider(start=minRef, end=maxRef, value=R, step=refStep,
                    title='$R_2$ (reflectivity of mirror 2)',
                    callback=callback,
                    #callback_policy = 'mouseup'
                    )
etaSlider = Slider(start=1.0, end=2.0, value=eta, step=0.01,
                   title = '$\eta$, refractive index of gas in etalon',
                   callback=callback)

callback.args['length'] = lengthSlider
callback.args['ref1'] = refl1Slider
callback.args['ref2'] = refl2Slider
callback.args['etaSl'] = etaSlider

div = Div(text="""
    <div> Free spectral range $\Delta\omega_{{FSR}}$ (in GHz): </div> <div id="FSR_info" class="valueDiv"></div>
"""#.format(deltaOMEGA_FSR*unitScale), width=400, height = 20)
          )

l = layout([
  [p1, p2],
  [widgetbox(refl1Slider, refl2Slider), widgetbox(div, lengthSlider, etaSlider)],
], sizing_mode='stretch_both')
output_file("airy.html", title="Fabry-Perot Cavity Intensity")
show(l)
