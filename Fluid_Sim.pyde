import math


N = 64
iter = 1
scale = 25

def constrain(x,mmin,mmax):
    return min(max(mmin,x),mmax)

def IX(x, y):
    x = constrain(x,0,N-1)
    y = constrain(y,0,N-1)
    return x+y*N


def settings():
  size(N*scale, N*scale);

class FluidSquare:
   # size: int
    #dt: float
    #diff: float
    #visc: float

    #s: np.array
    #density: np.array

    #Vx: np.array
    #Vy: np.array

    #Vx0: np.array
    #Vy0: np.array

    def __init__(self,  diffusion,  viscosity,  dt):

        self.size = N
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        self.s = [0] * (N*N)
        self.density = [0] * (N*N)

        self.Vx = [0] * (N*N)
        self.Vy = [0] * (N*N)

        self.Vx0 = [0] * (N*N)
        self.Vy0 = [0] * (N*N)

    def addDye(self, x,  y,  amount):
        index = IX(x, y)
        self.density[index] += amount

    def addVelocity(self, x,  y,  amountX, amountY):
        index = IX(x, y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY

    def step(self):
        visc = self.visc
        diff = self.diff
        dt = self.dt
        Vx = self.Vx
        Vy = self.Vy
        Vx0 = self.Vx0
        Vy0 = self.Vy0
        s = self.s
        density = self.density

        diffuse(1, Vx0, Vx, visc, dt)
        diffuse(2, Vy0, Vy, visc, dt)

        project(Vx0, Vy0, Vx, Vy)

        advect(1, Vx, Vx0, Vx0, Vy0, dt)
        advect(2, Vy, Vy0, Vx0, Vy0, dt)

        project(Vx, Vy, Vx0, Vy0)

        diffuse(0, s, density, diff, dt)
        advect(0, density, s, Vx, Vy, dt)

    def renderD(self):
        for i in range(0, N):
            for j in range(0, N):
                x = i * scale
                y = j * scale
                d = self.density[IX(i, j)]
                fill(d)
                noStroke()
                square(x,y,scale)
    def fadeD(self):
        for i in range(0,len(self.density)):
            d = self.density[i]
            self.density[i] = constrain(d-0.1,0,255)
                
def set_bnd(b, x):
    for i in range(1, N-1):
        x[IX(i, 0)] = -x[IX(i, 1)] if b == 2 else x[IX(i, 1)]
        x[IX(i, N-1)] = -x[IX(i, N-2)] if b == 2 else x[IX(i, N-2)]

    for j in range(1, N-1):
        x[IX(0, j)] = -x[IX(1, j)] if b == 1 else x[IX(1, j)]
        x[IX(N-1, j)] = -x[IX(N-2, j)] if b == 1 else x[IX(N-2, j)]

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)])
    x[IX(0, N-1)] = 0.5 * (x[IX(1, N-1)] + x[IX(0, N-2)])
    x[IX(N-1, 0)] = 0.5 * (x[IX(N-2, 0)] + x[IX(N-1, 1)])
    x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)])

def lin_solve(b, x, x0, a, c):
    cRecip = 1.0 / c
    for k in range(0, iter):
        for j in range(1, N-1):
            for i in range(1, N-1):
                x[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i+1, j)] + x[IX(i-1, j)] + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip

        set_bnd(b, x)


    
def project(velocX, velocY, p, div):
    for j in range(1, N-1):
        for i in range(1, N-1):
            div[IX(i, j)] = -0.5*(velocX[IX(i+1, j)] - velocX[IX(i-1, j)
                                                              ] + velocY[IX(i, j+1)] - velocY[IX(i, j-1)])/N
            p[IX(i, j)] = 0

    set_bnd(0, div)
    set_bnd(0, p)
    lin_solve(0, p, div, 1, 6)

    for j in range(1, N-1):
        for i in range(1, N-1):
            velocX[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]) * N
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]) * N

    set_bnd(1, velocX)
    set_bnd(2, velocY)
    
def diffuse(b, x, x0, diff, dt):
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a)
    
        
def advect(b, d, d0,  velocX, velocY,  dt):
    # i0, i1, j0, j1

    dtx = dt * (N - 2)
    dty = dt * (N - 2)

   # s0, s1, t0, t1
   # tmp1, tmp2, x, y

    Nfloat = N
    # ifloat, jfloat, kfloat
    # i, j

    for j, jfloat in zip(range(1, N-1), range(1, N-1)):
        for i, ifloat in zip(range(1, N-1), range(1, N-1)):
            tmp1 = dtx * velocX[IX(i, j)]
            tmp2 = dty * velocY[IX(i, j)]
            x = ifloat - tmp1
            y = jfloat - tmp2

            if(x < 0.5):
                x = 0.5
            if(x > Nfloat + 0.5):
                x = Nfloat + 0.5
            i0 = math.floor(x)
            i1 = i0 + 1.0
            if(y < 0.5):
                y = 0.5
            if(y > Nfloat + 0.5):
                y = Nfloat + 0.5
            j0 = math.floor(y)
            j1 = j0 + 1.0

            s1 = x - i0
            s0 = 1.0 - s1
            t1 = y - j0
            t0 = 1.0 - t1

            i0i = int(i0)
            i1i = int(i1)
            j0i = int(j0)
            j1i = int(j1)

            # review this line
            a = s0 * (t0 * d0[IX(i0i, j0i)]) + (t1 * d0[IX(i0i, j1i)])
            b = s1 * (t0 * d0[IX(i1i, j0i)]) + (t1 * d0[IX(i1i, j1i)])
            d[IX(i, j)] = a + b
        
    set_bnd(b, d)
    
    
fluid = FluidSquare(0,0,0.05)

def mouseDragged():
    fluid.addDye(mouseX/scale,mouseY/scale,500)
    amtX = mouseX - pmouseX
    amtY = mouseY - pmouseY
    fluid.addVelocity(mouseX/scale,mouseY/scale,amtX,amtY)

    
def draw():
    background(0)
    fluid.step()
    fluid.renderD()
    fluid.fadeD()
