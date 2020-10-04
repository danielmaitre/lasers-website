import numpy as np

# if js doesn't work, add
# <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjs/3.14.1/math.js"></script> 
# to the header of the html file

from bokeh.layouts import row, widgetbox, layout
from bokeh.models import CustomJS, Slider, RadioButtonGroup, Div
from bokeh.plotting import figure, output_file, show, ColumnDataSource

c = 3e8 # Speed of light, in m/s
L = 0.5 # cavity length, in m. Vary this for different results.
deltaOmega = 2*np.pi*c/(2*L) # mode spacing in angular frequency
offsetFrequency = 0 # vary this to add an overall shift to figure
E0 = 1

T = 2*L / c # round trip time
Ns = np.arange(0,100)
phis=np.zeros((100))
omegas = offsetFrequency + Ns * deltaOmega

maxMode=100
gran=3000 # granularity of plot, > 10000 takes too long in the browser to be interactive

#need to convert this to JS
def calcModelockedIntensity(N=100, modeLocking=True):
    Ns = np.arange(0,N) # generate array of integers from 0 to N
    omegas = offsetFrequency + Ns * deltaOmega # mode frequencies

    if modeLocking == True: # modelocking in operation, i.e. all phases identical
        phis = np.zeros((N)) # can arbitrarily choose the phase, set to zero for simplicity
    else: # otherwise, each mode has its own random phase
        phis = np.random.rand(N) * 2*np.pi

    
    
    ts = np.linspace(0,3,gran) # time, to be plotted on x-axis, in units of T

    E0 = 1 # maximum electric field amplitude

    Is = np.zeros((ts.size)) # will store the intensities at each time interval
    for i in range(ts.size): # for each instant in time, calculate laser intensity value
        tInSeconds = ts[i]*T
        Is[i] = np.absolute(E0*np.sum(np.exp(1j*(omegas*tInSeconds + phis))))**2 # sums all phase factors at time t due to the different modes

    return ts, Is#need to return max and mean intensity as well

x, y = calcModelockedIntensity()#change to take input from a slider, and add input from a toggle button or a tab
plot = figure(x_range=(0,3),
              title="Modelocking Demonstration",
              plot_width=400, plot_height=400)
plot.xaxis.axis_label = "time t, in units of cavity round trip time T"
plot.yaxis.axis_label = "Intensity in watts per sqaure metre"#could try adding latex labels



# write callback code and link needed variables
source = ColumnDataSource(data=dict(x=x, y=y, phis=np.append(phis,[0]*(gran-phis.size)),
                                    varbs=[T]+[offsetFrequency]+[deltaOmega]+[E0]+[0]*(gran-4)))
#print(x.size)

toggleSrce = ColumnDataSource(data=dict(phis=np.append(phis, [0]*(maxMode-phis.size))))

plot.line('x','y', source=source, line_width=1, line_alpha=1.0)

# this is a Javascript implementation of the calcModelockedIntensity function above
# it is called when the slider is changed
callback1 = CustomJS(args=dict(source=source, source2=toggleSrce), code = """ 
    // get the provided arrays from the python code
    var data = source.data;
    ts = data['x'];
    Is = data['y'];
    varbs = data['varbs'];

    var data2 = source2.data;
    phis = data2['phis'];
    
    //console.log("In callback, Is initial size, ts size, Is final size:");
    //console.log(Is.length);

    // store other variables from python code
    T = varbs[0];
    offsetFrequency = varbs[1];
    deltaOmega = varbs[2];
    E0 = varbs[3];

    var lockChoice = locking.active;
    var N = modes.value; // get current slider position

    if (lockChoice == 0){
        for (i=0; i<N; i++){
            phis[i] = 0.0; // experiment with matrices
        }
    } else {
        for (i=0; i<N; i++){
            phis[i] = math.random(2*math.pi);
        }
    }
    
    
    var omegas = new Array(N); // calculate the mode frequencies
    for (i=0; i<N; i++) {
        omegas[i] = offsetFrequency + i*deltaOmega;
    }

    IArray = new Array(0);
    for (j=0; j<ts.length; j++) { // calculate times to plot and the intensity values at those times using physics eqns
        ts[j] = j * 3/ts.length;

        var E=0;
        for (i=0; i<N; i++) {
            E = math.add(E, math.exp(math.complex(0,(omegas[i]*ts[j]*T + phis[i])))); // electric field value of mode i
        }
        E = math.multiply(E, E0);

        Is[j] = math.square(math.abs(E));
        IArray.push(Is[j]);
        //console.log(Is[j])
    }

    //console.log(ts.length);
    //console.log(Is.length);
    

    // output max intensity and mean intensity
    var maxDiv = document.getElementById('maxInfo');
    var meanDiv = document.getElementById('meanInfo');
    maxDiv.textContent = math.format(math.max(IArray),3);
    meanDiv.textContent = math.format(math.mean(IArray),3);

    source.change.emit(); // update the arrays
    source2.change.emit(); // update phis
""")

# this is a Javascript implementation of the if-else statement in the calcModelockedIntensity function above
# it is called when the radio button group is changed and recalculates phis
callback2 = CustomJS(args=dict(source1=source, source2=toggleSrce), code = """
    var data1 = source1.data;
    ts = data1['x'];
    Is = data1['y'];
    varbs = data1['varbs'];
    
    // store other variables from python code
    T = varbs[0];
    offsetFrequency = varbs[1];
    deltaOmega = varbs[2];
    E0 = varbs[3];
    
    var data2 = source2.data;
    phis = data2['phis'];
    
    var lockChoice = locking.active;
    var N = modes.value;

    if (lockChoice == 0){
        for (i=0; i<N; i++){
            phis[i] = 0.0; // experiment with matrices
        }
    } else {
        for (i=0; i<N; i++){
            phis[i] = math.random(2*math.pi);
        }
    }
    

    var omegas = new Array(N); // calculate the mode frequencies
    for (i=0; i<N; i++) {
        omegas[i] = offsetFrequency + i*deltaOmega;
    }

    IArray = new Array(0);
    for (j=0; j<ts.length; j++) { // calculate times to plot and the intensity values at those times using physics eqns
        ts[j] = j * 3/ts.length;

        var E=0;
        for (i=0; i<N; i++) {
            E = math.add(E, math.exp(math.complex(0,(omegas[i]*ts[j]*T + phis[i])))); // electric field value of mode i
        }
        E = math.multiply(E, E0);

        Is[j] = math.square(math.abs(E));
        IArray.push(Is[j]);
    }
    // output max intensity and mean intensity
    var maxDiv = document.getElementById('maxInfo');
    var meanDiv = document.getElementById('meanInfo');
    maxDiv.textContent = math.format(math.max(IArray),3);
    meanDiv.textContent = math.format(math.mean(IArray),3);

    source1.change.emit(); // update the arrays
    source2.change.emit(); // update phis
""")



# define widgets
modes_slider = Slider(start=2, end=maxMode, value=100, step=1, title = "Number of Modes, $N$", callback=callback1,
#                     callback_policy="mouseup"
                      )
lock_toggle = RadioButtonGroup(labels=["Modelocking", "No Modelocking"], active=0, callback=callback2)

callback1.args["modes"] = modes_slider
callback1.args["locking"] = lock_toggle
callback2.args["modes"] = modes_slider
callback2.args["locking"] = lock_toggle

maxI = np.amax(y)
meanI = np.mean(y)
print(meanI)

div1 = Div(text="""
    <div> Peak intensity of the current plot: </div> <div id="maxInfo" class="valueDiv"></div>
"""
           #.format(maxI), width=300, height=20 #<div id="initialMax">{0:3g}</div> removed [similar for mean below]
           )
div2 = Div(text="""
    <div> Mean intensity of the current plot: </div> <div id="meanInfo" class="valueDiv"></div>
"""
           #.format(meanI), width=300, height=20
           )

l = layout([[plot], [widgetbox(modes_slider, lock_toggle, div1, div2)]], sizing_mode='stretch_both')
output_file("modelocking.html", title="Modelocking")
show(l)







'''

if modeLocking ==True:
    addStr = ''  
else:
    addStr = 'non-'
            
maxI = np.amax(Is)
print('Maximum Intensity for this distribution is {0:.3g}'.format(maxI))
meanI = np.mean(Is)
print('Mean Intensity for this distribution is {0:.3g}'.format(meanI))
        
plt.figure()
plt.title('Intensity distribution produced over time by a {}modelocked laser\nwith {} modes and cavity length L = {} m'
    .format(addStr,N,L))
plt.xlabel('Time, in units of T (the round trip time = {0:.3g} s)'.format(T))
plt.ylabel('Intensity') # need units
plt.plot(ts, Is)
plt.show()
'''
