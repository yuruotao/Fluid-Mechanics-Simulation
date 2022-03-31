import os

import numpy as np


# CLASS CONFIGURATION 配置类
class Boundary:
    def __init__(self, boundary_type, boundary_value):
        self.DefineBoundary(boundary_type, boundary_value)

    def DefineBoundary(self, boundary_type, boundary_value):
        self.type = boundary_type
        self.value = boundary_value


# 设置边界条件，便于PDE的求解
# 种类有两种：dirichlet与neumann

class Space:
    def __init__(self):
        pass

    def CreateMesh(self, rowpts, colpts):
        self.rowpts = rowpts
        self.colpts = colpts
        self.u = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.v = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.p = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.p_c = np.zeros((self.rowpts, self.colpts))
        self.u_c = np.zeros((self.rowpts, self.colpts))
        self.v_c = np.zeros((self.rowpts, self.colpts))
        self.SetSourceTerm()

    def SetDeltas(self, breadth, length):
        self.dx = length / (self.colpts - 1)
        self.dy = breadth / (self.rowpts - 1)

    def SetInitialU(self, U):
        self.u = U * self.u

    def SetInitialV(self, V):
        self.v = V * self.v

    def SetInitialP(self, P):
        self.p = P * self.p

    def SetSourceTerm(self, S_x=0, S_y=0):
        self.S_x = S_x
        self.S_y = S_y


# 定义边界围成的区域


class Fluid:
    def __init__(self, rho, mu):
        self.SetFluidProperties(rho, mu)

    def SetFluidProperties(self, rho, mu):
        self.rho = rho
        self.mu = mu


# 定义区域内部的流体，其特征由密度和粘度表示


# 设置边界条件
# U,V为两轴，下面设置速度的边界条件

def SetUBoundary(space, left, right, top, bottom):
    if left.type == "D":
        space.u[:, 0] = left.value
    elif left.type == "N":
        space.u[:, 0] = -left.value * space.dx + space.u[:, 1]

    if right.type == "D":
        space.u[:, -1] = right.value
    elif right.type == "N":
        space.u[:, -1] = right.value * space.dx + space.u[:, -2]

    if top.type == "D":
        space.u[-1, :] = 2 * top.value - space.u[-2, :]
    elif top.type == "N":
        space.u[-1, :] = -top.value * space.dy + space.u[-2, :]

    if bottom.type == "D":
        space.u[0, :] = 2 * bottom.value - space.u[1, :]
    elif bottom.type == "N":
        space.u[0, :] = bottom.value * space.dy + space.u[1, :]


def SetVBoundary(space, left, right, top, bottom):
    if left.type == "D":
        space.v[:, 0] = 2 * left.value - space.v[:, 1]
    elif left.type == "N":
        space.v[:, 0] = -left.value * space.dx + space.v[:, 1]

    if right.type == "D":
        space.v[:, -1] = 2 * right.value - space.v[:, -2]
    elif right.type == "N":
        space.v[:, -1] = right.value * space.dx + space.v[:, -2]

    if top.type == "D":
        space.v[-1, :] = top.value
    elif top.type == "N":
        space.v[-1, :] = -top.value * space.dy + space.v[-2, :]

    if bottom.type == "D":
        space.v[0, :] = bottom.value
    elif bottom.type == "N":
        space.v[0, :] = bottom.value * space.dy + space.v[1, :]


# 压力边界条件
def SetPBoundary(space, left, right, top, bottom):
    if left.type == "D":
        space.p[:, 0] = left.value
    elif left.type == "N":
        space.p[:, 0] = -left.value * space.dx + space.p[:, 1]

    if right.type == "D":
        space.p[1, -1] = right.value
    elif right.type == "N":
        space.p[:, -1] = right.value * space.dx + space.p[:, -2]

    if top.type == "D":
        space.p[-1, :] = top.value
    elif top.type == "N":
        space.p[-1, :] = -top.value * space.dy + space.p[-2, :]

    if bottom.type == "D":
        space.p[0, :] = bottom.value
    elif bottom.type == "N":
        space.p[0, :] = bottom.value * space.dy + space.p[1, :]


# 设置函数
def SetTimeStep(CFL, space, fluid):
    with np.errstate(divide='ignore'):
        dt = CFL / np.sum([np.amax(space.u) / space.dx, np.amax(space.v) / space.dy])
    # Escape condition if dt is infinity due to zero velocity initially
    if np.isinf(dt):
        dt = CFL * (space.dx + space.dy)
    space.dt = dt


# 时间间隔函数


def GetStarredVelocities(space, fluid):
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u = space.u.astype(float, copy=False)
    v = space.v.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    S_x = float(space.S_x)
    S_y = float(space.S_y)
    rho = float(fluid.rho)
    mu = float(fluid.mu)

    u_star = u.copy()
    v_star = v.copy()

    u1_y = (u[2:, 1:cols + 1] - u[0:rows, 1:cols + 1]) / (2 * dy)
    u1_x = (u[1:rows + 1, 2:] - u[1:rows + 1, 0:cols]) / (2 * dx)
    u2_y = (u[2:, 1:cols + 1] - 2 * u[1:rows + 1, 1:cols + 1] + u[0:rows, 1:cols + 1]) / (dy ** 2)
    u2_x = (u[1:rows + 1, 2:] - 2 * u[1:rows + 1, 1:cols + 1] + u[1:rows + 1, 0:cols]) / (dx ** 2)
    v_face = (v[1:rows + 1, 1:cols + 1] + v[1:rows + 1, 0:cols] + v[2:, 1:cols + 1] + v[2:, 0:cols]) / 4
    u_star[1:rows + 1, 1:cols + 1] = u[1:rows + 1, 1:cols + 1] - dt * (
            u[1:rows + 1, 1:cols + 1] * u1_x + v_face * u1_y) + (dt * (mu / rho) * (u2_x + u2_y)) + (dt * S_x)

    v1_y = (v[2:, 1:cols + 1] - v[0:rows, 1:cols + 1]) / (2 * dy)
    v1_x = (v[1:rows + 1, 2:] - v[1:rows + 1, 0:cols]) / (2 * dx)
    v2_y = (v[2:, 1:cols + 1] - 2 * v[1:rows + 1, 1:cols + 1] + v[0:rows, 1:cols + 1]) / (dy ** 2)
    v2_x = (v[1:rows + 1, 2:] - 2 * v[1:rows + 1, 1:cols + 1] + v[1:rows + 1, 0:cols]) / (dx ** 2)
    u_face = (u[1:rows + 1, 1:cols + 1] + u[1:rows + 1, 2:] + u[0:rows, 1:cols + 1] + u[0:rows, 2:]) / 4
    v_star[1:rows + 1, 1:cols + 1] = v[1:rows + 1, 1:cols + 1] - dt * (
            u_face * v1_x + v[1:rows + 1, 1:cols + 1] * v1_y) + (dt * (mu / rho) * (v2_x + v2_y)) + (dt * S_y)

    space.u_star = u_star.copy()
    space.v_star = v_star.copy()


# 速度


def SolvePressurePoisson(space, fluid, left, right, top, bottom):
    # Save object attributes as local variable with explicit typing for improved readability
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u_star = space.u_star.astype(float, copy=False)
    v_star = space.v_star.astype(float, copy=False)
    p = space.p.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    rho = float(fluid.rho)
    factor = 1 / (2 / dx ** 2 + 2 / dy ** 2)

    error = 1
    tol = 1e-3

    ustar1_x = (u_star[1:rows + 1, 2:] - u_star[1:rows + 1, 0:cols]) / (2 * dx)
    vstar1_y = (v_star[2:, 1:cols + 1] - v_star[0:rows, 1:cols + 1]) / (2 * dy)

    i = 0
    while error > tol:
        i += 1
        p_old = p.astype(float, copy=True)
        p2_xy = (p_old[2:, 1:cols + 1] + p_old[0:rows, 1:cols + 1]) / dy ** 2 + (
                p_old[1:rows + 1, 2:] + p_old[1:rows + 1, 0:cols]) / dx ** 2
        p[1:rows + 1, 1:cols + 1] = p2_xy * factor - (rho * factor / dt) * (ustar1_x + vstar1_y)
        error = np.amax(abs(p - p_old))
        # Apply Boundary Conditions
        SetPBoundary(space, left, right, top, bottom)

        if i > 500:
            tol *= 10


# 泊松压力方程


def SolveMomentumEquation(space, fluid):
    # Save object attributes as local variable with explicit typing for improved readability
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u_star = space.u_star.astype(float)
    v_star = space.v_star.astype(float)
    p = space.p.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    rho = float(fluid.rho)
    u = space.u.astype(float, copy=False)
    v = space.v.astype(float, copy=False)

    p1_x = (p[1:rows + 1, 2:] - p[1:rows + 1, 0:cols]) / (2 * dx)
    u[1:rows + 1, 1:cols + 1] = u_star[1:rows + 1, 1:cols + 1] - (dt / rho) * p1_x

    p1_y = (p[2:, 1:cols + 1] - p[0:rows, 1:cols + 1]) / (2 * dy)
    v[1:rows + 1, 1:cols + 1] = v_star[1:rows + 1, 1:cols + 1] - (dt / rho) * p1_y


# 动量方程


def SetCentrePUV(space):
    space.p_c = space.p[1:-1, 1:-1]
    space.u_c = space.u[1:-1, 1:-1]
    space.v_c = space.v[1:-1, 1:-1]


# 保存边界内的速度与压力


def MakeResultDirectory(wipe=False):
    cwdir = os.getcwd()
    dir_path = os.path.join(cwdir, "Result")
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path, exist_ok=True)
    else:
        if wipe:
            os.chdir(dir_path)
            filelist = os.listdir()
            for file in filelist:
                os.remove(file)

    os.chdir(cwdir)


# 将结果保存到Result里


def WriteToFile(space, iteration, interval):
    if iteration % interval == 0:
        dir_path = os.path.join(os.getcwd(), "Result")
        filename = "PUV{0}.txt".format(iteration)
        path = os.path.join(dir_path, filename)
        with open(path, "w") as f:
            for i in range(space.rowpts):
                for j in range(space.colpts):
                    f.write("{}\t{}\t{}\n".format(space.p_c[i, j], space.u_c[i, j], space.v_c[i, j]))
# 在迭代间隔中将变量值进行保存
