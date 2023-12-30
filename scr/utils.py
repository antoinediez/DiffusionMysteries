import os
import sys
import cv2
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import EllipseCollection


def make_video(number_of_frames=None,simu_name="simu",video_name='funny_video',rate=40,prefix="t_",every=1):
    img_array = []
    current_directory = os.getcwd()
    frame_directory = current_directory+"/"+simu_name+"/frames"
    if number_of_frames is None:
        number_of_frames = len([name for name in os.listdir(frame_directory) if os.path.isfile(os.path.join(frame_directory, name))])
    for count in range(number_of_frames):
        if count%every==0:
            filename = frame_directory+"/"+prefix+str(count)+'.png'
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width,height)
            img_array.append(img)
    out = cv2.VideoWriter(simu_name+"/"+video_name+'.avi',cv2.VideoWriter_fourcc(*'DIVX'), rate,size)
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()
    
def run_hardspheres(simu,T,t_iter,plot_every,simu_name,record_traj=False,plot_traj=False,trajcolor='k',trajwidth=1.8):
    percent = 0
    T0 = simu.t
    
    simu.update_plot()
    simu.fig.savefig(simu_name + "/frames/" + f"{t_iter}.png")
    t_iter += 1
    
    plot_time = 0.0
    
    if record_traj or plot_traj:
        traj = {'x' : [simu.pos[0,0].item()],
                'y' : [simu.pos[0,1].item()]}
        traj_plot, = simu.ax.plot(traj['x'],traj['y'],color='k',linewidth=trajwidth)
    else:
        traj = None
        traj_plot = None
            
    for _ in simu:
        plot_time += simu.dt
        if record_traj or plot_traj:
            traj['x'].append(simu.pos[0,0].item())
            traj['y'].append(simu.pos[0,1].item())
        if (simu.t - T0)/(T-T0) >= percent/100:
            sys.stdout.write('\r' + "Progress:" + str(percent) + "%")
            sys.stdout.flush()
            percent += 1
        if plot_time>=plot_every:
            simu.update_plot()
            if plot_traj:
                traj_plot.set_xdata(traj['x'])
                traj_plot.set_ydata(traj['y'])
            simu.fig.savefig(simu_name + "/frames/" + f"{t_iter}.png")
            t_iter += 1
            plot_time = 0.0
        if simu.t>T:
            break
        
    extra = {"traj" : traj,
            "traj_plot" : traj_plot}
            
    return t_iter, extra

def run_hardspheres_tp(simu1,simu2,parent_fig,T,t_iter,plot_every,simu_name):
    percent = 0
    T0 = simu1.t
    
    simu1.update_plot()
    simu2.update_plot()
    parent_fig.savefig(simu_name + "/frames/" + f"{t_iter}.png")
    t_iter += 1
    
    plot_time = 0.0
    
    # if record_traj or plot_traj:
    #     traj = {'x' : [simu.pos[0,0].item()],
    #             'y' : [simu.pos[0,1].item()]}
    #     traj_plot, = simu.ax.plot(traj['x'],traj['y'],color='k',linewidth=trajwidth)
    # else:
    #     traj = None
    #     traj_plot = None
            
    for _ in simu2:
        simu1.__next__()
        plot_time += simu2.dt
        # if record_traj or plot_traj:
        #     traj['x'].append(simu.pos[0,0].item())
        #     traj['y'].append(simu.pos[0,1].item())
        if (simu2.t - T0)/(T-T0) >= percent/100:
            sys.stdout.write('\r' + "Progress:" + str(percent) + "%")
            sys.stdout.flush()
            percent += 1
        if plot_time>=plot_every:
            simu1.update_plot()
            simu2.update_plot()
            # if plot_traj:
            #     traj_plot.set_xdata(traj['x'])
            #     traj_plot.set_ydata(traj['y'])
            parent_fig.savefig(simu_name + "/frames/" + f"{t_iter}.png")
            t_iter += 1
            plot_time = 0.0
        if simu2.t>T:
            break
        
    # extra = {"traj" : traj,
    #         "traj_plot" : traj_plot}
    extra = None
            
    return t_iter, extra
    
    
# def scatter_particles_onehighlight(simu, time, show=True,
#                       save=False, path='simu'):
#     """Scatter plot with the radii of the particles.

#     Args:
#         simu (Particles): A model.
#         time (list): List of times to plot.
#         show (bool, optional): Show the plot. Default is True.
#         save (bool, optional): Save the plots. Default is False.
#         path (str, optional): The plots will be saved in the directory
#             ``'./path/frames'``.

#     """
#     # Scatter plots at times specified by the list 'time' for the
#     # mechanism 'mechanism' (iterator inherited from the class 'particle')
#     t = 0
#     tmax = len(time) - 1
#     percent = 0
#     real_time = [0.]
#     x0 = [simu.pos[0,0].item()]
#     y0 = [simu.pos[0,1].item()]
#     if save:
#         current_directory = os.getcwd()
#         final_directory = os.path.join(current_directory, path, r'frames')
#         if not os.path.exists(final_directory):
#             os.makedirs(final_directory)
#     for it, info in enumerate(simu):
#         real_time.append(real_time[-1] + simu.dt)
#         if real_time[-1] / time[tmax] >= percent / 100:
#             sys.stdout.write('\r' + "Progress:" + str(percent) + "%")
#             sys.stdout.flush()
#             percent += 1
#         if abs(real_time[-1] >= time[t]):
#             if simu.d == 2:
#                 f = plt.figure(t, figsize=(6, 6))
#                 ax = f.add_subplot(111)
#                 x = simu.pos[:, 0].cpu()
#                 y = simu.pos[:, 1].cpu()
#                 x0.append(simu.pos[0,0].item())
#                 y0.append(simu.pos[0,1].item())
#                 traj = ax.plot(x0,y0,'k',linewidth=3)
#                 size = 2 * simu.R.cpu()
#                 offsets = list(zip(x, y))
#                 fcolors = ["tab:blue"] * simu.N
#                 fcolors[0] = "tab:orange"
#                 op = 0.5*np.ones(simu.N)
#                 op[0] = 1.0
#                 ax.add_collection(EllipseCollection(
#                     widths=size, heights=size, angles=0, units='xy',facecolor=fcolors,
#                     edgecolor='k', offsets=offsets, transOffset=ax.transData,alpha=op))
#                 ax.axis('equal')  # set aspect ratio to equal
#                 ax.axis([0, simu.L[0].cpu(), 0, simu.L[1].cpu()])
#                 #                 ax.set_aspect(1)
#                 #                 plt.axes().set_aspect('equal', 'box')
#                 ax.set_title(simu.name + '\n Parameters: ' + simu.parameters
#                              + '\n Time=' + str(round(real_time[-1], 1)),
#                              fontsize=10)
#             elif simu.d == 3:
#                 f = plt.figure(t, figsize=(10, 10))
#                 ax = f.add_subplot(111, projection='3d')
#                 X = simu.pos[:, 0].cpu()
#                 Y = simu.pos[:, 1].cpu()
#                 Z = simu.pos[:, 2].cpu()
#                 ax.scatter(X, Y, Z, s=0.1)
#                 ax.set_xlim(0, simu.L[0].cpu())
#                 ax.set_ylim(0, simu.L[1].cpu())
#                 ax.set_zlim(0, simu.L[2].cpu())
#                 ax.set_title(
#                     simu.name + '\n Parameters: ' + simu.parameters
#                     + '\n Time=' + str(round(simu.iteration * simu.dt, 1)),
#                     fontsize=10)
#             if save:
#                 f.savefig(f"{final_directory}/" + str(t) + ".png")
#             if show:
#                 # plt.show()
#                 plt.pause(0.000001)
#             else:
#                 plt.close()
#             t += 1
#         if t > tmax:
#             break

