''' Disentangled Airy plotting from specHoles_finesse.py'''

import numpy as np
from bokeh.layouts import layout, row, column, widgetbox
from bokeh.models import CustomJS, Slider, Dropdown, Div
from bokeh.plotting import figure, output_file, show, ColumnDataSource

''' constants '''
pi = np.pi
c = 3e8
ln2 = np.log(2)


''' Laser cavity variables (all changeable via a slider) '''
L = 0.2 # length of cavity in m
alpha = 1.0 # absorption and scattering losses in m^-1
R = 0.99 # reflectivity of mirrors
offset = 0; # frequency offset of zeroth mode from zero detuning

inputPower = 10.0
maxInputPower = 10.0

T = 1200 + 273 # discharge T in kelvin (check whether function of atom or not)

''' Atom variables (dropdown menu to pick the atom) '''
Gamma = 2*pi*5e6  #*450e6 # linewidth of each mode - needed to make it this small to avoid overlapping
                            #Lorentzian effects which distorted the results ( a proper correction in the
                            #code will allow for larger. more realistic values)
lambd = 488e-9 # laser wavelength
M_A = 40 # relative atomic mass of Ar

''' Calculated variables (will need to recalculate when one of the slider variables changes) '''
omega_12 = 2*pi*c / lambd # laser transition frequency
deltaOMEGA_FSR = pi*c/L # Free spectral range mode spacing in cavity (reducing L increases spacing)
delta_c = np.log(1/(R**2)) + 2*alpha*L # fractional round trip loss (increasing reflectivities reduces this
                                            #without affecting mode spacing)
deltaOMEGA_I = 7.16e-7 * omega_12 * (T / M_A)**(0.5) # linewidth of the inhomogeneously (Doppler) broadened Gaussian

lzFWHM = Gamma

def gaussian(det):
    return 2/(deltaOMEGA_I) * (ln2/pi)**(0.5) * np.exp(-4*(det**2) * ln2 / (deltaOMEGA_I**2))
def lorentzian(det):
    return (lzFWHM/(2*pi)) / (det**2 + (lzFWHM**2)/4)

stickHeight = -0.1 # mainly arbitrary, how high to make the mode sticks (just enough to be noticeable)
stickBottom = -0.1*maxInputPower

''' plotting variables '''
gran = 1000 # how detailed to make plot (reduce to increase speed)
unitScale = 1/(2*pi*1e9) # Ghz
plotWidth = 1.5 # how much of gaussian to show

''' work out the maximum number of modes possible in the cavity '''
maxLength = 3.0 # slider maximum
min_deltaOMEGA_FSR = pi*c/maxLength
maxModes = int(2*deltaOMEGA_I*plotWidth / min_deltaOMEGA_FSR) + 1


''' ---------------- CREATE PLOT ARRAYS (will need to ve implemented in JS callbacks) ------------- '''

''' update the following when offset slider changes '''
moddedOffset1 = offset % deltaOMEGA_FSR # calculate closest mode peak to right of resonance
moddedOffset2 = deltaOMEGA_FSR - moddedOffset1 # calculate closest mode peak to left of resonance

''' following arrays only need to be changed if atom changes (i.e. only for dropdown box) '''

detunings = np.linspace(-deltaOMEGA_I*plotWidth, deltaOMEGA_I*plotWidth, gran)
detGHz = detunings * unitScale # scale to plotting units
Ls = np.zeros((gran)) # Lineshape due to inhomogeneous Doppler broadening
for i in range(gran):
    Ls[i] = gaussian(detunings[i])
maxL = np.max(Ls)
Ls *= inputPower/maxL # scale peak to inputPower


''' ****** update the following when anything changes  ****** '''
n = 0
Lz1=np.zeros((gran)) # stores the current Lorentzian
lzTot=np.zeros((gran)) # contains the hole-burned lineshape
bar_positions = np.array([]) # frequency location for the output power spikes (also serves as modePositions)
P_outs = np.array([])
print(bar_positions)
print(P_outs)

stickHeights = np.array([])
cont1=True # set to False when the Doppler lineshape drops below delta_c or we reach end of plot (right side)
cont2=True # set to False when the Doppler lineshape drops below delta_c or we reach end of plot (left side)
while (cont1==True | cont2==True ): # while there are still Lorentzians to subtract from the total lineshape
    ''' find where the current Lorentzians are centred relative to the gaussian '''
    centredFreq = n*deltaOMEGA_FSR 
    loc1 = centredFreq + moddedOffset1 # right side
    loc2 = -centredFreq - moddedOffset2 # left side
    
    if ( cont1 == True & (inputPower*gaussian(loc1)/maxL > delta_c)
         & ((loc1) < (deltaOMEGA_I*plotWidth + lzFWHM*plotWidth))):  # try to make second limit as small as possible for speed
        for i in range(gran):
            Lz1[i] = lorentzian(detunings[i]-loc1) # Lorentzian
        Lz1 /= np.max(Lz1)
        for i in range(gran):
            gap = Ls[i] - lzTot[i] - delta_c
            if gap > 0:
                Lz1[i] = Lz1[i]*gap
            else:
                Lz1[i] = 0
        lzTot += Lz1

        if abs(loc1) != abs(loc2): # mirror the Lorentzian to the other side as long as it doesnt't already
                                    # coincide with a mode
            if ((inputPower*gaussian(-loc1)/maxL > delta_c) & ((-loc1) > (-deltaOMEGA_I*plotWidth - lzFWHM*plotWidth))):
                for i in range(gran):
                    detuning = detunings[i]
                    Lz1[i] = lorentzian(detunings[i]+loc1) # Lorentzian
                Lz1 = Lz1/np.max(Lz1)
                for i in range(gran):
                    gap = Ls[i] - lzTot[i] - delta_c
                    if gap > 0:
                        Lz1[i] = Lz1[i]*gap
                    else:
                        Lz1[i] = 0
                lzTot += Lz1

        # for each mode, add to power plot
        P = (inputPower*gaussian(loc1)/maxL)/delta_c - 1
        bar_positions = np.append(bar_positions, loc1*unitScale)
        stickHeights = np.append(stickHeights, stickHeight)
        P_outs = np.append(P_outs, P) # in arbitrary units, we know P is proportional to g_0/g_th - 1
    else:
        cont1=False
        bar_positions = np.append(bar_positions, loc1*unitScale)
        stickHeights = np.append(stickHeights, stickHeight) # add sticks even if not present in power spectrum
        P_outs = np.append(P_outs, 0) # so doesn't appear in power graph

    if ( cont2 == True & (inputPower*gaussian(loc2)/maxL > delta_c)
         & ((loc2) > (-deltaOMEGA_I*plotWidth - lzFWHM*plotWidth))):  # try to make second limit as small as possible for speed
        for i in range(gran):
            Lz1[i] = lorentzian(detunings[i]-loc2) # Lorentzian
        Lz1 = Lz1/np.max(Lz1)
        for i in range(gran):
            gap = Ls[i] - lzTot[i] - delta_c
            if gap > 0:
                Lz1[i] = Lz1[i]*gap
            else:
                Lz1[i] = 0
        lzTot += Lz1

        if abs(loc1) != abs(loc2): # mirror
            if ((inputPower*gaussian(-loc2)/maxL > delta_c) & ((-loc2) < (deltaOMEGA_I*plotWidth + lzFWHM*plotWidth))):
                for i in range(gran):
                    Lz1[i] = lorentzian(detunings[i]+loc2) # Lorentzian
                Lz1 = Lz1/np.max(Lz1)
                for i in range(gran):
                    gap = Ls[i] - lzTot[i] - delta_c
                    if gap > 0:
                        Lz1[i] = Lz1[i]*gap
                    else:
                        Lz1[i] = 0
                lzTot += Lz1

        P = (inputPower*gaussian(loc2)/maxL)/delta_c - 1
        bar_positions = np.append(bar_positions, loc2*unitScale)
        stickHeights = np.append(stickHeights, stickHeight)
        P_outs = np.append(P_outs, P) # in arbitrary units, we know P is proportional to g_0/g_th - 1
    else:
        cont2=False
        bar_positions = np.append(bar_positions, loc2*unitScale)
        stickHeights = np.append(stickHeights, stickHeight) # add sticks even if not present in power spectrum
        P_outs = np.append(P_outs, 0) # so doesn't appear in power graph
  
    n+=1

cont1=True
cont2=True
while (cont1==True | cont2==True ): # add additional sticks to left graph up to edge
    centredFreq = n*deltaOMEGA_FSR 
    loc1 = centredFreq + moddedOffset1 # right side
    loc2 = -centredFreq - moddedOffset2 # left side

    if ((cont1==True) & ( loc1 < deltaOMEGA_I*plotWidth)):
        bar_positions = np.append(bar_positions, loc1*unitScale)
        stickHeights = np.append(stickHeights, stickHeight) # add sticks even if not present in power spectrum
        P_outs = np.append(P_outs, 0) # so doesn't appear in power graph
    else:
        cont1 = False

    if ((cont2==True) & ( loc2 > -deltaOMEGA_I*plotWidth)):
        bar_positions = np.append(bar_positions, loc2*unitScale)
        stickHeights = np.append(stickHeights, stickHeight) # add sticks even if not present in power spectrum
        P_outs = np.append(P_outs, 0) # so doesn't appear in power graph
    else:
        cont2 = False

    n+=1

for i in range(P_outs.size):
    if P_outs[i] < 0:
        P_outs[i] = 0
    
''' ***************************************************** '''

''' ----------------------------------------------------------------------------------------------- '''


''' callback data sources '''
atomAttr = ColumnDataSource(data=dict(Gamma=[Gamma], dOmegaI=[deltaOMEGA_I], refl=[R],
                                      plotWidth=[plotWidth], maxL=[maxL], stickH=[stickHeight],
                                      sBottom=[stickBottom], scale=[unitScale], maxModes=[maxModes])) # never updated unless atom changed #, M_A=[M_A], omega12=[omega_12]
plot1Arrs = ColumnDataSource(data=dict(x=detunings, y=(Ls - lzTot), Ls=(Ls/inputPower), det=detGHz))
plotDeltaC = ColumnDataSource(data=dict(x=[detGHz[0], detGHz[gran-1]], y=[delta_c]*2)) # updated when delta_c is
plot2Arrs = ColumnDataSource(data=dict(x=np.append(bar_positions,[0]*(maxModes - bar_positions.size)),
                                       y=np.append(P_outs,[0]*(maxModes - P_outs.size)),
                                       sh=np.append(stickHeights,[stickBottom]*(maxModes - stickHeights.size))))

''' ************** create plots **************** '''
p1 = figure(x_range=(-deltaOMEGA_I*plotWidth*unitScale,deltaOMEGA_I*plotWidth*unitScale),
            y_range=(stickBottom,maxInputPower),
            title='Hole-burned gain distribution in frequency',
            plot_width=400, plot_height=400)
p1.xaxis.axis_label = "detuning from laser transition frequency, in GHz"
p1.yaxis.axis_label = 'round trip gain (arb units)'
p1.line('det', 'y', source=plot1Arrs, line_width=1, line_alpha=1.0, color='blue', legend="Round trip gain")
p1.line('x', 'y', source=plotDeltaC, line_width=1, line_alpha=1.0, color='red', legend="Fractional round trip loss")
p1.vbar(x='x', width=None, top='sh', bottom=stickBottom, source=plot2Arrs, color='black', line_dash='dashed',
        legend="Cavity Mode Positions")
p1.legend.location = "top_left"
p1.legend.background_fill_alpha = 0.5


p2 = figure(x_range=(-deltaOMEGA_I*plotWidth*unitScale,deltaOMEGA_I*plotWidth*unitScale), title='Output Power at different frequencies',
            plot_width=400, plot_height=400)
p2.xaxis.axis_label = "detuning from laser transition frequency, in GHz"
p2.yaxis.axis_label = "laser output power (arb units)"
p2.vbar(x='x', width=None, top='y', bottom=0.0, source=plot2Arrs, color='black')

''' ******************************************* '''


''' ---------------------------------- CALLBACKS ------------------------------------------- '''
lengthCb = CustomJS(args=dict(srce1=atomAttr, srce4=plot1Arrs, srce5=plotDeltaC,
                              srce6=plot2Arrs), code = """
    var c = 3e8; // speed of light
    /****************** GET ALL DATA **********************************/
    var data1 = srce1.data;
    Gamma = data1['Gamma'][0]; // change [0] indices to value of dropdown box if added
    deltaOMEGA_I = data1['dOmegaI'][0];
    R = data1['refl'][0];
    plotWidth = data1['plotWidth'][0];
    maxL = data1['maxL'][0];
    sHeight = data1['stickH'][0];
    stickBottom = data1['sBottom'][0];
    unitScale = data1['scale'][0];
    maxModes = data1['maxModes'][0];

    var data4 = srce4.data;
    detunings = data4['x'];
    gainOut = data4['y'];
    Ls = data4['Ls'];
    gran = detunings.length;

    var data5 = srce5.data;
    flatH = data5['y'];

    var data6 = srce6.data;
    barPositions = data6['x'];
    P_outs = data6['y'];
    stickHeights = data6['sh'];

    /* widgets */
    var L = length.value;
    var alpha = alph.value;
    var offset = off.value;
    var inputPower = pow.value;
    /****************************************************************/
    
    /* recalculate free spectral range spacing */
    var deltaOMEGA_FSR = math.pi*c/L;

    /* update modded offsets */
    var moddedOffset1 = offset % deltaOMEGA_FSR;
    var moddedOffset2 = deltaOMEGA_FSR - moddedOffset1;

    /* recalculate delta_c */
    var delta_c = math.log(1/(R*R)) + 2*alpha*L;
    
    var lzFWHM = Gamma;
    
    /* update all plots */
    var n = 0;
    var lz = new Array(gran);
    var totArray = new Array(gran);
    for (var i = 0; i < gran; i++) {
        totArray[i] = 0;
    }
    var newBarPositions = new Array(0);
    var newPs = new Array(0);
    var newShs = new Array(0);
    var cont1 = true; var cont2 = true;
    while ((cont1==true) || (cont2==true)){
        var centredFreq = n*deltaOMEGA_FSR;
        var loc1 = centredFreq + moddedOffset1;
        var loc2 = -centredFreq - moddedOffset2;
        
        if ((cont1==true) && (inputPower*gaussian(loc1, deltaOMEGA_I)/maxL > delta_c) && ( loc1 < (deltaOMEGA_I*plotWidth + lzFWHM*plotWidth))) {
            for (var i=0; i < gran; i++) {
                lz[i] = lorentzian(detunings[i]-loc1, lzFWHM);
            }
            var maxes = math.max(lz);
            lz = math.divide(lz, maxes);
            for (var i=0; i < gran; i++) {
                gap = inputPower*Ls[i] - totArray[i] - delta_c;
                if (gap > 0)
                    lz[i] = lz[i]*(gap);
                else
                    lz[i] = 0;
            }
            totArray = math.add(totArray, lz);

            if (math.abs(loc1) != math.abs(loc2)) {
                for (var i=0; i < gran; i++) {
                    lz[i] = lorentzian(detunings[i]+loc1, lzFWHM);
                }
                maxes = math.max(lz);
                lz = math.divide(lz, maxes);
                for (var i=0; i < gran; i++) {
                    gap = inputPower*Ls[i] - totArray[i] - delta_c;
                    if (gap > 0)
                        lz[i] = lz[i]*(gap);
                    else
                        lz[i] = 0;
                }
                totArray = math.add(totArray, lz); 
            }
            P = (inputPower*gaussian(loc1, deltaOMEGA_I)/maxL)/delta_c - 1;
            newBarPositions.push(loc1*unitScale);
            newShs.push(sHeight);
            newPs.push(P);
        } else {
            cont1 = false;
            newBarPositions.push(loc1*unitScale);
            newShs.push(sHeight);
            newPs.push(0);
        }

        if ((cont2==true) && (inputPower*gaussian(loc2, deltaOMEGA_I)/maxL > delta_c) && ( loc2 > (-deltaOMEGA_I*plotWidth - lzFWHM*plotWidth))) {
            for (var i=0; i < gran; i++) {
                lz[i] = lorentzian(detunings[i]-loc2, lzFWHM);
            }
            maxes = math.max(lz);
            lz = math.divide(lz, maxes);
            for (var i=0; i < gran; i++) {
                gap = inputPower*Ls[i] - totArray[i] - delta_c;
                if (gap > 0)
                    lz[i] = lz[i]*(gap);
                else
                    lz[i] = 0;
            }
            totArray = math.add(totArray, lz);

            if (math.abs(loc1) != math.abs(loc2)) {
                for (var i=0; i < gran; i++) {
                    lz[i] = lorentzian(detunings[i]+loc2, lzFWHM);
                }
                maxes = math.max(lz);
                lz = math.divide(lz, maxes);
                for (var i=0; i < gran; i++) {
                    gap = inputPower*Ls[i] - totArray[i] - delta_c;
                    if (gap > 0)
                        lz[i] = lz[i]*(gap);
                    else
                        lz[i] = 0;
                }
                totArray = math.add(totArray, lz); 
            }
            P = (inputPower*gaussian(loc2, deltaOMEGA_I)/maxL)/delta_c - 1;
            newBarPositions.push(loc2*unitScale);
            newShs.push(sHeight);
            newPs.push(P);
        } else {
            cont2 = false;
            newBarPositions.push(loc2*unitScale);
            newShs.push(sHeight);
            newPs.push(0);
        }
        
        n++;       
    }

    /* add additional sticks to left graph up to edge */
    cont1 = true; cont2 = true;
    while ((cont1==true) || (cont2==true)){
        var centredFreq = n*deltaOMEGA_FSR;
        var loc1 = centredFreq + moddedOffset1;
        var loc2 = -centredFreq - moddedOffset2;

        if ((cont1==true) && ( loc1 < deltaOMEGA_I*plotWidth)) {
            newBarPositions.push(loc1*unitScale);
            newShs.push(sHeight);
            newPs.push(0);
        } else { cont1 = false; }

        if ((cont2==true) && ( loc2 > -deltaOMEGA_I*plotWidth)) {
            newBarPositions.push(loc2*unitScale);
            newShs.push(sHeight);
            newPs.push(0);
        } else { cont2 = false; }

        n++;
    }

    flatH[0] = delta_c; flatH[1] = delta_c;
    for (var i=0; i<gran; i++) {
        gainOut[i] = inputPower*Ls[i] - totArray[i]; 
    }
    for (var i=0; i < newBarPositions.length; i++) {
        barPositions[i] = newBarPositions[i];
        P_outs[i] = newPs[i];
        stickHeights[i] = newShs[i];
    }
    for (var i = newBarPositions.length; i < maxModes; i++) {
        P_outs[i] = 0; // stop other bars from being plotted
        stickHeights[i] = stickBottom;
    }

    /* update all changed data */
    srce4.change.emit(); // update gainOut
    srce5.change.emit(); // update flatH
    srce6.change.emit(); // update bar_positions, stickHeights and P_outs

""")



''' ---------------------------------------------------------------------------------------- '''

# need to define the gaussian, lorentzian and maxAndIndex functions used in the Custom JS in the html.
# In the bat file, need to escape angles and closing round-brackets i.e. < becomes ^<, > becomes ^>, ) becomes ^)
gaussian = """
    function gaussian(det, fwhm) {
        return 2/(fwhm) * math.pow((math.log(2)/math.pi),0.5) * math.exp((-4*math.log(2)*math.square(det)) / math.square(fwhm));
    }
"""

lorentzian = """
    function lorentzian(det, lw) {
        return (lw/(2*math.pi)) / (math.square(det) + math.square(lw)/4);
    }
"""


''' ------------------------------------- WIDGETS ------------------------------------------ '''
powerSlider = Slider(start=0.5, end=maxInputPower, value=inputPower,
                      step=0.5, title="input power (arb. units)",
                      callback=lengthCb,
                      )
lengthSlider = Slider(start=0.01, end=maxLength, value=L, step=0.01,
                      title='$L$ (cavity length in m)',
                      callback=lengthCb,
                      )
alphaSlider = Slider(start = 0.005, end=2.000, value=alpha, step=0.005,
                     title='$\\alpha$ (absorption and scattering losses in cavity in m$^{-1}$)',
                     callback=lengthCb,
                     )
offsetSlider = Slider(start=-deltaOMEGA_I*plotWidth, end=deltaOMEGA_I*plotWidth, value=offset,
                      step=2*deltaOMEGA_I*plotWidth/(gran), title="modes offset in s$^{-1}$",
                      callback=lengthCb,
                      )


lengthCb.args["pow"] = powerSlider
lengthCb.args["length"] = lengthSlider
lengthCb.args["alph"] = alphaSlider
lengthCb.args["off"] = offsetSlider

''' ---------------------------------------------------------------------------------------- '''


''' create layout and show '''
l = layout([
  [p1, p2],
  [widgetbox(powerSlider, lengthSlider, alphaSlider, offsetSlider)], 
], sizing_mode='stretch_both')
output_file("holeBurning.html", title="Spectral Hole Burning")
show(l)
