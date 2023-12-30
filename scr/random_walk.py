import numpy as np

class RandomWalk:
    
    def __init__(self,x=np.array([[5.,5.]]),rate=1.0,dt=0.01,falling_rate=42.0,compute_y=True,bound=None,bc="periodic"):
        self.x = x
        self.N = self.x.shape[0]
        self.dim = self.x.ndim
        self.dt = dt
        self.rate = rate
        self.bound = bound
        self.bc = bc
        self.compute_y = compute_y
        if self.dim==1:
            if self.compute_y:
                count = self.count(self.x)
                self.y = np.zeros(self.N)
                for x in range(int(np.min(self.x)),int(np.max(self.x)+1)):
                    count = (self.x==x).sum()
                    self.y[self.x==x] = np.arange(0,count,1)
                self.falling_rate = falling_rate
                self.to_jump = np.zeros(self.N)
                self.after_jump = np.zeros(self.N)
        else:
            self.compute_y = False ## Not Implemented...         
                        
    def update(self):
        who_jumps = np.random.rand(self.N) > np.exp(-self.dt*self.rate)
        if self.dim==2:
            jump = (np.random.rand(who_jumps.sum(),1)<0.5)
            jump = np.concatenate((jump,1.0-jump),axis=1)
            jump *= (2.0*(np.random.rand(who_jumps.sum(),1)<0.5) - 1.0)
            self.x[who_jumps,:] += jump
        elif self.dim==1:
            if self.compute_y:
                for i in range(self.N):
                    if who_jumps[i] and self.to_jump[i]==0.0:
                        self.to_jump[i] = 1.0
                    if self.to_jump[i]>0.0:
                        self.to_jump[i] += 1.0
                    if self.to_jump[i]==5.0:
                        new_x = self.x[i] + 2.0*(np.random.rand()<0.5) - 1.0
                        self.to_jump[i] = 0.0
                        self.after_jump[i] = 1.0
                        self.y[i] = np.max(self.y[self.x==new_x])+1.0 if (self.x==new_x).sum()>0 else 0.0
                        self.x[i] = new_x
                    if self.after_jump[i]>0.0:
                        self.after_jump[i] += 1.0
                    if self.after_jump[i]==5.0:
                        self.after_jump[i]=0.0
                        
                for x in range(int(np.min(self.x)),int(np.max(self.x)+1)):
                    count = (self.x==x).sum()
                    y_target = np.zeros(count)
                    y_target[np.argsort(self.y[self.x==x])] = np.arange(0,count,1)
                    if self.falling_rate>0:
                        self.y[self.x==x] = np.maximum(y_target,self.y[self.x==x] - ((y_target-self.y[self.x==x])<0.0)*self.falling_rate*self.dt)
                    else:
                        self.y[self.x==x] = y_target
            else:
                self.x[who_jumps] += (2.0*(np.random.rand(who_jumps.sum())<0.5) - 1.0)
        if not self.bound is None:
            if self.bc == "periodic":
                self.x = self.x % self.bound
            elif self.bc == "reflecting":
                self.x[self.x>=self.bound] = self.bound - 1
                self.x[self.x<0] = 0
            
    def count(self,K):
        if self.dim==1:
            return np.array([(self.x==k).sum() for k in K])
        else:
            raise NotImplementedError()
        
    def histogram(self,K1,K2):
        if self.dim>1:
            raise NotImplementedError()
        bins = np.arange(start=K1-0.5,stop=K2+0.5,step=1.0)
        return np.histogram(self.x,bins=bins)
    
    def add_particles(self,x):
        self.x = np.concatenate((self.x,x),axis=0)
        self.N += len(x)
        if self.dim==1:
            if self.compute_y:
                count = self.count(self.x)
                self.y = np.zeros(self.N)
                for x in range(int(np.min(self.x)),int(np.max(self.x)+1)):
                    count = (self.x==x).sum()
                    self.y[self.x==x] = np.arange(0,count,1)
                self.to_jump = np.zeros(self.N)
                self.after_jump = np.zeros(self.N)