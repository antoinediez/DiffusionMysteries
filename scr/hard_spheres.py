import math
import numpy as np 
import torch
from pykeops.torch import LazyTensor
import sisyphe
from matplotlib import pyplot as plt
from matplotlib.collections import EllipseCollection

from sisyphe.particles import KineticParticles
from sisyphe.kernels import squared_distances

class HardSpheres(KineticParticles):
    
    def __init__(self, pos, vel, radius, box_size, boundary_conditions, dt0, fig=None):
        
        super().__init__(pos=pos, vel=vel, 
                         interaction_radius=radius, 
                         box_size=box_size, 
                         boundary_conditions=boundary_conditions)
        
        self.dt0 = dt0
        self.dt = dt0
        self.t = 0.0
        self.radius = radius
        self.R = radius * torch.ones(self.N,device=self.pos.device,dtype=self.pos.dtype) if isinstance(radius,float) else radius
        self.mass = self.R ** 2
        self.name = "Hard Spheres"
        self.parameters = 'N=' + str(self.N)
        self.fig = plt.figure(figsize=(8,8)) if fig is None else fig
        self.ax = self.fig.add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.hs_plot = self.__plot_hardspheres()        
    
    @property
    def momentum(self):
        return self.vel.sum(0)/self.N
    
    @property
    def kinetic_energy(self):
        return 0.5*(self.vel ** 2).sum()/self.N
        
    def good_config(self):
        
        sq_dist = squared_distances(
        self.pos,
        self.pos,
        self.L,
        boundary_conditions = self.bc)
        
        R_i = LazyTensor(self.R.reshape((self.pos.shape[0], 1))[:, None])
        R_j = LazyTensor(self.R.reshape((self.pos.shape[0], 1))[None, :])
        R_ij = R_i + R_j
        
        K = (R_ij ** 2 - sq_dist).relu()
        # bad_guys = (K.sum(1) - (2*self.radius)**2).reshape(self.N)
        bad_guys = K.sum(1).reshape(self.N) - (2*self.R)**2
        good = (bad_guys.sum() == 0.0)
        
        if good:
            return good, None
        else:
            col_index = (-K).argKmin(2,1)
            col_index = col_index[:,1]
            col_index[bad_guys==0] = -1.0
            return False, col_index
        
        
    def solve_col(self,col_index):
        pos_col1 = self.pos[col_index>=0,:]
        vel_col1 = self.vel[col_index>=0,:]
        pos_col2 = self.pos[col_index[col_index>=0],:]
        vel_col2 = self.vel[col_index[col_index>=0],:]
        r_col1 = self.R[col_index>=0]
        r_col2 = self.R[col_index[col_index>=0]]
        m_col1 = self.mass[col_index>=0]
        m_col2 = self.mass[col_index[col_index>=0]]
        
        dx = pos_col1 - pos_col2
        dv = vel_col1 - vel_col2
        r = r_col1 + r_col2
        
        dt_pre = ((dx*dv).sum(1) + ((dx*dv).sum(1)**2 - (dv*dv).sum(1) * ((dx*dx).sum(1) - r**2)).sqrt())/((dv*dv).sum(1))
        dt_pre = dt_pre.reshape((dt_pre.size()[0],1))
        pos_col1_pre = pos_col1 - dt_pre * vel_col1
        # pos_col2_pre = pos_col2 - dt_pre * vel_col2
        
        nu12 = dx / torch.norm(dx,dim=1).reshape(dx.size()[0],1)
        s = (nu12 * dv).sum(1).reshape(dx.size()[0],1) * nu12
        vel1_post = vel_col1 - s * (2*m_col2/(m_col1+m_col2)).reshape(dx.size()[0],1)
        vel2_post = vel_col1 + vel_col2 - vel1_post
        
        kinet12_pre = m_col1 * (vel_col1 ** 2).sum(1) + m_col2 * (vel_col2 ** 2).sum(1)
        kinet12_post = m_col1 * (vel1_post ** 2).sum(1) + m_col2 * (vel2_post ** 2).sum(1)
        vel1_post *= kinet12_pre.sqrt().reshape(dv.size()[0],1)/kinet12_post.sqrt().reshape(dv.size()[0],1)
        
        kinet0 = 0.5 * (m_col1.reshape(dv.size()[0],1) * (vel_col1 ** 2)).sum()
        kinet_post =  0.5 * (m_col1.reshape(dv.size()[0],1) * (vel1_post ** 2)).sum()
        
        self.pos[col_index>=0,:] = pos_col1_pre + dt_pre*vel1_post
        # self.pos[col_index[col_index>=0],:] = pos_col2_pre + dt_pre*vel2_post
        self.vel[col_index>=0,:] = vel1_post * kinet0.sqrt()/kinet_post.sqrt()  # scale the velocity to preserve kinetic energy without error
        # self.vel[col_index[col_index>=0],:] = vel2_post
        # Note: it is a loop so update only one particle. 
    
    def check_boundary_walls(self):
        # center = torch.tensor([[0.5,0.5]],device=self.pos.device,dtype=self.pos.dtype)
        # out = torch.max((self.pos - center).abs(),1) > (0.5-self.R)
        for i in range(self.d):
            if self.bc[i] == 1:
                inf0 = self.pos[:, i] < self.R
                supL = self.pos[:, i] > self.L[i] - self.R
                self.pos[inf0, i] = -self.pos[inf0, i] + 2*self.R[inf0]
                # self.pos[supL, i] = 2 * self.L[i] - self.pos[supL, i] + self.R[supL]
                self.pos[supL, i] -=  2*(self.pos[supL, i] - (self.L[i] - self.R[supL]))

                self.vel[inf0, i] = -self.vel[inf0, i]
                self.vel[supL, i] = -self.vel[supL, i]
            elif self.bc[i] == 0:
                self.pos[:, i] = torch.remainder(self.pos[:, i], self.L[i])
        
        
    def update(self):
        pos0 = self.pos.detach().clone()
        vel0 = self.vel.detach().clone()
        kinet0 = self.kinetic_energy
        self.dt = 2*self.dt0
        is_good = False
        while (not(is_good) and (self.dt > 0.000001)):
            self.dt /= 2
            self.pos = pos0 + vel0*self.dt
            self.check_boundary_walls()
            is_good, col_index = self.good_config()
            if ~is_good:
                self.solve_col(col_index)
            is_good, col_index = self.good_config()
        self.t += self.dt
        
    def __plot_hardspheres(self):
        size = 2 * self.R.cpu()
        x = self.pos[:, 0].cpu()
        y = self.pos[:, 1].cpu()
        offsets = list(zip(x, y))
        fcolors = ["tab:blue"] * self.N
        op = 0.7*np.ones(self.N)
        hs_plot = self.ax.add_collection(EllipseCollection(
            widths=size, heights=size, angles=0, units='xy',facecolor=fcolors,
            edgecolor='k', offsets=offsets, transOffset=self.ax.transData,alpha=op))
        self.ax.axis('equal')  # set aspect ratio to equal
        self.ax.axis([0, self.L[0].cpu(), 0, self.L[1].cpu()])
        return hs_plot
    
    def update_plot(self):
        x = self.pos[:, 0].cpu()
        y = self.pos[:, 1].cpu()
        offsets = list(zip(x, y))
        self.hs_plot.set_offsets(offsets)