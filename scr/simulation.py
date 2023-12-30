import numpy as np
import matplotlib as mpl
from matplotlib.collections import EllipseCollection
from matplotlib.collections import PolyCollection


class Simulation:
    
    def __init__(self,rw,lat,part_size=0.8,clock_rate=1.0,N_nbgh_plot=1,N_traj_plot=1,traj_linewidth=None,traj_alpha=0.8):
        self.rw = rw
        self.lat = lat
        if self.lat.dim != self.rw.dim:
            raise ValueError("The dimensions of the lattice and the random walk must be equal.")
        self.dim = self.lat.dim
        self.N_nbgh_plot = N_nbgh_plot
        self.part_size = part_size
        self.particles, self.extra = self.__plot_particles()
        self.histo = self.__plot_histo()
        if self.rw.compute_y:
            self.histo.set_alpha(0.0)
        self.clock_rate = clock_rate
        self.traj = [self.rw.x.copy()]
        self.N_traj_plot = N_traj_plot
        self.trajs_plot = self.__plot_traj(linewidth=traj_linewidth,alpha=traj_alpha) if self.N_traj_plot>0 else None
        
    
    def update(self,update_nbgh=True,update_traj=True):
        self.rw.update()
        
        if self.dim==2 or self.rw.compute_y: 
            offsets = self.__compute_pos_offsets()
            self.particles.set_offsets(offsets)
        
        if self.dim==1:
            if self.rw.compute_y:
                offsets = self.__compute_pos_offsets()
                self.particles.set_offsets(offsets)
                fcolors = ["tab:red" if b else "tab:orange" for b in (self.rw.to_jump+self.rw.after_jump)>0]
                self.particles.set_facecolor(fcolors)
            else:
                count,bins = self.rw.histogram(K1=self.lat.L/2.-self.lat.zoom/2,K2=self.lat.L/2.+self.lat.zoom/2)
                self.histo.set(data=count)
            
        if update_nbgh and self.N_nbgh_plot>0:
            self.update_neighbours()
        
        if update_traj:
            if self.N_traj_plot>0:
                self.traj.append(self.rw.x.copy())
                
            for n in range(self.N_traj_plot):
                x = [pos[n,0] for pos in self.traj]
                y = [pos[n,1] for pos in self.traj]
                self.trajs_plot[n].set_xdata(x)
                self.trajs_plot[n].set_ydata(y)
            
        
        self.lat.set_time(self.lat.time+self.clock_rate*self.rw.dt)
    
    def update_neighbours(self):
        verts = self.__compute_neighbours()
        self.extra["nbgh"].set_verts(verts)
       
    def __plot_histo(self,fill=True,facecolor="tab:orange"):
        if self.dim>1:
            return None
        count,bins = self.rw.histogram(K1=self.lat.L/2.-self.lat.zoom/2,K2=self.lat.L/2.+self.lat.zoom/2)
        return self.lat.ax.stairs(values=count,edges=bins,fill=fill,facecolor=facecolor)
        
    
    def __plot_traj(self,linewidth=None,alpha=0.8):
        if linewidth is None:
            linewidth = 2*self.lat.x_lines.get_linewidth() 
        trajs = []
        for n in range(self.N_traj_plot):
            x = [pos[n,0] for pos in self.traj]
            y = [pos[n,1] for pos in self.traj]
            traj_plot, = self.lat.ax.plot(x,y,'k',linewidth=linewidth,alpha=alpha)
            trajs.append(traj_plot)
        return trajs
    
    
    def __plot_particles(self,fcolors="tab:orange"):
        if self.dim==1 and not self.rw.compute_y:
            return None,None
        offsets = self.__compute_pos_offsets()
        particles = self.lat.ax.add_collection(
                EllipseCollection(
                    widths=self.part_size, heights=self.part_size, angles=0, units='xy',facecolor=fcolors,
                    edgecolor='k', offsets=offsets, transOffset=self.lat.ax.transData
                )
            )
        extra = {}
        
        if self.dim==2:            
            if self.N_nbgh_plot>0:
                verts = self.__compute_neighbours()
                nbgh = self.lat.ax.add_collection(
                    PolyCollection(verts,alpha=0.5,zorder=-1.0)
                )
                extra["nbgh"] = nbgh
            
        
        return particles, extra
                
        
    def __compute_pos_offsets(self):
        if self.dim==2:
            x = self.rw.x[:, 0]
            y = self.rw.x[:, 1]
        elif self.dim==1:
            x = self.rw.x
            y = (self.rw.y+0.5)*self.part_size
        return list(zip(x, y))
        
    def __compute_neighbours(self):
        verts = []
        for n in range(self.N_nbgh_plot):
            nbgh_west = np.array(
                [[self.rw.x[n,0]-1.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]-0.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]-0.5,self.rw.x[n,1]+0.5],
                    [self.rw.x[n,0]-1.5,self.rw.x[n,1]+0.5],
                ]
            )
            verts.append(nbgh_west)
            
            nbgh_east = np.array(
                [[self.rw.x[n,0]+1.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]+0.5],
                    [self.rw.x[n,0]+1.5,self.rw.x[n,1]+0.5],
                ]
            )
            verts.append(nbgh_east)
            
            nbgh_north = np.array(
                [[self.rw.x[n,0]-0.5,self.rw.x[n,1]+1.5],
                    [self.rw.x[n,0]-0.5,self.rw.x[n,1]+0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]+0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]+1.5],
                ]
            )
            
            verts.append(nbgh_north)
            
            nbgh_south = np.array(
                [[self.rw.x[n,0]-0.5,self.rw.x[n,1]-1.5],
                    [self.rw.x[n,0]-0.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]-0.5],
                    [self.rw.x[n,0]+0.5,self.rw.x[n,1]-1.5],
                ]
            )
            
            verts.append(nbgh_south)
        return verts