from matplotlib import pyplot as plt
import numpy as np

class Lattice:
    
    def __init__(self,dim=2,L=10,clock_size=None,fig=None):
        if clock_size is None:
            clock_size = 0.25*L
        self.fig = plt.figure(figsize=(8,8)) if fig is None else fig
        if clock_size>0:
            gs = self.fig.add_gridspec(
                2,2,
                width_ratios=[1,clock_size/L],
                height_ratios=[clock_size/L,1],
                wspace = 0.0, hspace=0.0
            )
            self.ax = self.fig.add_subplot(gs[1,0])
            self.clock = self.fig.add_subplot(gs[0,1])
            self.clock.axis('equal')
            self.clock.axis("off")
            self.axtitle = self.fig.add_subplot(gs[0,0])
            self.axtitle.axis('equal')
            self.axtitle.axis("off")
            self.axtitle.text = self.axtitle.text(0.0,0.0,"Title",fontsize=16,ha="center",va="center")
            self.time = 1.0
            self.clock,self.hour,self.minute = self.__plot_clock()
        else:
            self.ax = self.fig.add_subplot(111)
            self.clock = None
            self.title = None
        # self.ax.axis('equal')
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.dim = dim
        self.L = L if L%2==0 else L+1
        self.zoom = L
        self.x_line_max = self.L+0.5 if self.dim==2 else 1.0 
        self.x_lines, self.y_lines = self.__plot_lattice()
        self.set_height()
        self.set_zoom(self.zoom)
        
        
    def __plot_lattice(self):
        mesh = np.arange(start=0.5,stop=self.L+0.5,step=1)
        x_lines = self.ax.vlines(mesh,-0.5,self.x_line_max,'k')
        y_lines = self.ax.hlines(mesh,-0.5,self.L+0.5,'k') if self.dim==2 else None
        return x_lines, y_lines
    
    def set_zoom(self,K):
        self.zoom = K if K%2==0 else K+1
        self.ax.set_xlim(self.L/2.-K/2.-0.45,self.L/2.+K/2.+0.45)
        if self.dim==2:
            self.ax.set_ylim(self.L/2.-K/2.-0.45,self.L/2+K/2.+0.45)
    
    def set_height(self,H=None):
        if self.dim==2:
            pass 
        else:
            if H is None:
                H = self.L/2
            self.ax.set_ylim(0.0,H)
            segs = self.x_lines.get_segments()
            for seg in segs:
                seg[1,1] = self.x_line_max
            self.set_lines_param(segments=segs)
                
    def set_lines_param(self,**params):
        self.x_lines.set(**params)
        if not self.y_lines is None:
            self.y_lines.set(**params)
            
    def __plot_clock(self):
        th = np.linspace(0,2*np.pi,1000)
        clock = self.clock.plot(0.5+0.5*np.cos(th),0.5+0.5*np.sin(th),'k')
        angle_hour = np.pi/2 - 2*np.pi*self.time/12
        angle_minute = np.pi/2 - 2*np.pi*(self.time%12)
        hour, = self.clock.plot([0.5,0.5+0.3*np.cos(angle_hour)],[0.5,0.5+0.3*np.sin(angle_hour)],'k',linewidth=3.0)
        minute, = self.clock.plot([0.5,0.5+0.4*np.cos(angle_minute)],[0.5,0.5+0.4*np.sin(angle_minute)],'k',linewidth=3.0)
        return clock,hour,minute
    
    def set_time(self,time):
        self.time = time 
        angle_hour = np.pi/2 - 2*np.pi*self.time/12
        angle_minute = np.pi/2 - 2*np.pi*(self.time%12)
        self.hour.set_xdata([0.5,0.5+0.3*np.cos(angle_hour)])
        self.hour.set_ydata([0.5,0.5+0.3*np.sin(angle_hour)])
        self.minute.set_xdata([0.5,0.5+0.4*np.cos(angle_minute)])
        self.minute.set_ydata([0.5,0.5+0.4*np.sin(angle_minute)])
        
    def set_title(self,title):
        self.axtitle.text.set(text=title)
        