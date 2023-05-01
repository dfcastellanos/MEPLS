import json
from io import StringIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
import argparse

def access_by_path(dic, path):   
    
    """
    
    Allows to access a nested dictionary using path to an entry, as if the nested
    dictionary was a file system.
    
    Parameters
    ----------
    - dic: nested dictionary
    - path: string with a path to a value within the dictionary (/key1/key2/key3)

    Returns
    -------
    The object located in the input path
    
    """
    
    keys = path.strip('/').split('/')
    for key in keys:
        dic = dic[key]

    return dic


def load_data(filename, path_to_plastic, path_to_driving):

    """
    
    Load data from the input JSON file into a pandas DataFrame combining the
    platic and driving event histories, located in the input paths within the
    JSON file structure.
    
    Parameters
    ----------
    - filename: input JSON file
    - path_to_plastic: path to the plastic events history within the JSON file
                       (e.g. /Data/histories/plastic_events)    
    - path_to_driving: path to the driving events history within the JSON file
        
    Returns
    -------
    The object located in the input path
    
    """

    # load the JSON file as a nested dictionary    
    sim_data = json.load( open(filename) )
    
    parameters = sim_data['Parameters']
    Nx = int(parameters['Nx'])
    Ny = int(parameters['Ny'])
    
    # in the JSON file, there are strings containing CSV. We can create a pandas DataFrames with 
    # those strings as if they were input file using StringIO. We use the simulation history index
    # as the DataFrame index
    plastic_events = access_by_path(sim_data, path_to_plastic)
    plastic_df = pd.read_csv( StringIO(plastic_events) ).set_index('index')
   
    driving_events = access_by_path(sim_data, path_to_driving)
    driving_df = pd.read_csv( StringIO(driving_events) ).set_index('index')
    
    # merge both dataframes. Except the index, the DataFrames don't contain each other's columns,
    # so the merged DataFrame will have many NaNs. Since the DataFrame containts stress and strain
    # increments, and we will calculate accumulated magniudes, we can just replace the NaNs
    # by zeros 
    df = plastic_df.merge(driving_df, how='outer', left_index=True, right_index=True)
    df.loc[:,:] = df.loc[:,:].fillna(0.)

    # compute the local von Mises plastic strain increments
    df['dstrain'] = np.sqrt((2.)*(0.5*np.power(df.eigenstrain_00.values-df.eigenstrain_11.values,2.)+2*np.power(df.eigenstrain_01.values,2.)))
   
    # compute the accumulated magnitudes
    df['vm_plastic_strain'] = np.cumsum(df.dstrain.values)/(Nx*Ny)

    if 'dext_stress' in driving_df.columns:
        df['ext_stress'] = np.cumsum(df.dext_stress.values)
    if 'dtotal_strain' in driving_df.columns:
        df['total_strain'] = np.cumsum(df.dtotal_strain.values)
    if 'dload' in driving_df.columns:
        df['load'] = np.cumsum(df.dload.values)
    if 'dtime' in driving_df.columns:
        df['time'] = np.cumsum(df.dtime.values)

    # make sure that the element nubmers are integers, since they will be used for index arrays
    df['element'] = df['element'].astype(int)

    return df, parameters


def animate_simulation(X, Y, Z, elements, Nx, Ny, x_label=None, y_label=None, z_label=None, fps=2, duration=25, p0=(0.,0.), title=None):

    """
    
    Create an animation using the input data. It shows a left plot with 
    the Y(X) cruve, and a right plot with a spatial map created using values 
    from Z and the locations given by the elements. It also shows a color bar 
    linked to the spatial map color scale.
    
    Parameters
    ----------
    - X: array with the X data
    - Y: array with the Y data
    - Z: array with the Z data, corresponding to local changes of some field
    - elements: array with the element numbers in which the local changes of Z occured
    - Nx: number of mesoscale elements forming the system in the x-direction
    - Ny: number of mesoscale elements forming the system in the y-direction
    - x_label: x-label for the Y(X) plot
    - y_label: y-label for the Y(X) plot
    - z_label: label under the colorbar for the 2D plot of Z
    - fps: frames-per-second in the animation
    - duration: duration of the animation in seconds
    - p0: a point to be added at the beginning of the (X,Y) data series
    - title: set the title of the figure

    Returns
    -------
    The animation object (needed, since it will play as long as a refernce to it exists)
    
    Interactive controls
    --------------------
    - press spacebar: pause/resume the animation
    - press i: swicth between incremental mode (displays activity in a certain 
               window of the x axis magnitude) and integrated mode (show accumulated activity)
    - press ,: reduce the length of the window
    - press .: increase the length of the window
    - press up arrow: if mode is incremental, move the window one step forward; 
                      if mode is integrated, play animation in forward direction
    - press down arrow: if mode is incremental, move the window one step backward;
                        if mode is integrated, play animation in backward direction
    - press p: print the simulations parameters. It needs an entry /Parameters in the 
               input JSON file.    
        
    """
    
    # make sure that the input are numpy arrays, since are extremely faster than pandas series when indexing them
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    elements = np.array(elements)
    
    # add the initial point to the curve
    if p0 is not None:
        X = np.insert(X, 0, p0[0])
        Y = np.insert(Y, 0, p0[1])
        # a dummy initial point to the others, otherwise their history index would be shited by 1 with respsect to X and Y
        # after adding the initial point p0 to them
        Z = np.insert(Z, 0, 0.)
        elements = np.insert(elements, 0, 0)
  
    # setup the figure: the left axis with the Y(X) curve, the right curve with a 2D spatial map  for Z
    fig, axes = plt.subplots(1,2,figsize=(9,4), facecolor='white')

    if title is not None: 
        fig.suptitle(title)
               
    axes[0].set_xlabel(x_label, fontsize=12)
    axes[0].set_ylabel(y_label, fontsize=12)
    axes[0].set_xlim([min(X)*0.9,max(X)*1.02])
    axes[0].set_ylim([min(Y),max(Y)*1.02])
    axes[0].grid(True)
    
    # show the full curve in grey
    line2, = axes[0].plot(X, Y, linewidth=0.5, color='grey', alpha=0.5)
    
    # create the line object that we will update in the animation, showing the curve in red
    line, = axes[0].plot(X, Y, linewidth=0.5, color='red')
    
    # create the line object that we will update in the animation, showing the 2D map of Z
    im = axes[1].imshow(np.zeros([Nx,Ny]), cmap='afmhot', animated=True, interpolation='gaussian')

    plt.sca(axes[1])
    plt.tick_params(top=False, right=False)
    
    # make a color bar with the values of the 2D map
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cax.set_xlabel(z_label)
    cb = fig.colorbar(im, cax=cax)

    plt.tight_layout()



    # define how the plots change with the animation frame number
    def animate(i):
                
        anim.current_i = i

        if anim.incremental: 
            # to show only the data over an interval of frames, centered in the current frame
            i0 = i - 10 * anim.h
            i1 = i + 10 * anim.h
            if i1 >= len(elements):
                i1 = len(elements)-1
        else: 
            # to show the data from the beginning till the current frame
            i0 = 0
            i1 = i
        
        # update the curve
        line.set_xdata(X[i0:i1])
        line.set_ydata(Y[i0:i1])
           
        # using the elements numbers, recreate the history of Z up to index i1
        map_deformation = np.zeros(Nx*Ny)    
        for i in range(i0,i1-1):
            j = elements[i]
            map_deformation[j] += Z[i]  
            
        # so far we were using MEPLS convention of 1D vectors, but now we need to 
        # transform the map into 2D
        map_deformation = map_deformation.reshape(Nx,Ny)
        

        # update the 2D plot and adjust the color scale limits
        im.set_data(map_deformation)
        vmax = np.max(map_deformation)
        vmin = np.min(map_deformation)
        im.set_clim(vmin, vmax)
        
        return line, im 
        
    
    
    # define the controls of the figure's interactive behavior
    def on_press(event):
        
        # click space bar -> pause / resume the animation
        if event.key.isspace():
            if anim.running:
                anim.event_source.stop()
            else:
                anim.event_source.start()
            anim.running ^= True

        # click down arrow -> backward time evolution
        if event.key == 'down':
            anim.direction = -1
            if not anim.running:
                i = next(anim.frame_seq)
                line, im = animate(i)
                fig.canvas.draw()
                
        # click up arrow -> forward time evolution
        if event.key == 'up':       
            anim.direction = +1
            if not anim.running:
                i = next(anim.frame_seq)
                line, im = animate(i)
                fig.canvas.draw()

        # click 'p' -> print to the screen the simulation parameters
        if event.key == 'p':
            print('Simulation parameters:\n', '\n'.join(['{} = {}'.format(k,v) for k,v in parameters.items()]))

        # click 'i' -> toggle incremental mode
        if event.key == 'i':
            anim.incremental ^= True

        # click ',' -> reduce the lenght of the x-axis window in incremental mode
        if event.key == '.':
            anim.h += 1

        # click '.' -> increase the lenght of the x-axis window in incremental mode
        if event.key == ',':
            anim.h -= 1    
    
        # when a press event is detected, refresh the plots to show the changes
        line, im = animate(anim.current_i)
        fig.canvas.draw()

        return

    # use the key press events defined above in the figure we created
    fig.canvas.mpl_connect('key_press_event', on_press)
    
    
    
    # define how we move from one animation frame to the next.
    # We use a generator instead of an array of frames so we can modify
    # its state on real time, and change the animation diretion (forward or backward)
    def update_time():  
        
        # animation frames are a subset of the data array index
        nframes = int(fps*duration)
        frames = np.linspace(0, len(X), nframes).astype(int)    
        
        i=0
        # each step of the animation corresponds to the next or the previous frame
        # depending on the animation direction (forward or backward)
        while i < len(frames):
            yield frames[i]
            i += anim.direction
            
        return
            
            
            
    # create the actual animation with the figure and functions defined above
    anim = animation.FuncAnimation(fig, animate, update_time, interval=1, repeat_delay=1e3, blit=0)
    
    # we can add new members to the existing animation object, that will be used by the functions above
    anim.running = True
    anim.direction = +1
    anim.incremental = False
    anim.h = 5
    anim.current_i = 0
    
    return anim



if __name__ == '__main__':    

    """
    Example of use:
    python ../animate.py --data out.json /Data/athermal/plastic_events /Data/athermal/driving_events 
                         --vars total_strain ext_stress dstrain 
                         --labels '$\varepsilon_{\rm xy}$' '$\Sigma_{\rm xy}$' r'$\varepsilon_{\rm vm}(\vec{r})$'
    """

    controls = """
                Interactive controls
                --------------------
                - press spacebar: pause/resume the animation
                - press i: swicth between incremental mode (displays activity in a certain 
                           window of the x axis magnitude) and integrated mode (show accumulated activity)
                - press ,: reduce the length of the window
                - press .: increase the length of the window
                - press up arrow: if mode is incremental, move the window one step forward; 
                                  if mode is integrated, play animation in forward direction
                - press down arrow: if mode is incremental, move the window one step backward;
                                    if mode is integrated, play animation in backward direction
                - press p: print the simulations parameters. It needs an entry /Parameters in the 
                           input JSON file.    
                """
    
    # parse the command line arguments      
    parser = argparse.ArgumentParser(description='Make an animation with the results of the input JSON file.', usage=controls)
    parser.add_argument('data', type=str, nargs=3, help='Path to the JSON output file and to the plastic and driving events with. Example -d out.json Data/athermal/plastic_events Data/athermal/driving_events')
    parser.add_argument('vars', type=str, nargs=3, help='Names of the x, y and z variables. One plot show the y(x) curve and the other a 2D image with z.')
    parser.add_argument('-t', '--title', type=str, help='Set the title of the figure')
    parser.add_argument('-l', '--labels', type=str, nargs='+', help='Axies labels for x, y and z. Latex mode can be used between $ symbols.')
    parser.add_argument('-r', '--rescale', type=float, nargs='+', help='Factors to rescale the values of x, y and z respectively.')
    parser.add_argument('-o', '--output', type=str, help='Filename to save the animation as a gif')
    parser.add_argument('-i', '--dpi', type=int, help='DPI for the gif animation')

    args = vars(parser.parse_args())
        
    filename = args['data'][0]
    path_plastic_events = args['data'][1]
    path_driving_events = args['data'][2]
    x_name = args['vars'][0]
    y_name = args['vars'][1]
    z_name = args['vars'][2]
    
    try: x_label = args['labels'][0].replace('//', '/')
    except: x_label = x_name
    try: y_label = args['labels'][1].replace('//', '/')
    except: y_label = y_name
    try: z_label = args['labels'][2].replace('//', '/')
    except: z_label = z_name
    
    try: ax = args['rescale'][0]
    except: ax = 1
    try: ay = args['rescale'][1]
    except: ay = 1
    try: az = args['rescale'][2]
    except: az = 1    

        
    # load and prepare the data for the animation    
    df, parameters = load_data(filename, path_plastic_events, path_driving_events)
    
    X = df[x_name] * ax
    Y = df[y_name] * ay
    Z = df[z_name] * az
    
    elements = df['element']
    Nx = int(parameters['Nx'])
    Ny = int(parameters['Ny'])
        
        
    # make the animation
    anim = animate_simulation(X, Y, Z, elements, Nx, Ny, x_label, y_label, z_label, title=args['title'])    
    # without this call, the animation won't show if the Python session is not interactive
    
    
    # to save the animation as a gif
    if args['output'] is not None:
        # anim.save(args['output'], writer='imagemagick', dpi=args['dpi'] or 70)
        # to save the animation as a mp4 video
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=10, bitrate=1800)
        anim.save('out.mp4', writer=writer, dpi=200)
    else:
        plt.show()
