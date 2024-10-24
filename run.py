from ctypes import CDLL, POINTER, c_int, c_double
import taichi as ti
import numpy as np
from numpy.ctypeslib import ndpointer

nd_double_array_dim1 = ndpointer(dtype=np.double, ndim=1, flags='C')
nd_double_array_dim2 = ndpointer(dtype=np.double, ndim=2, flags='C')
nd_int_array_dim1 = ndpointer(dtype=np.int32, ndim=1, flags='C')
nd_int_array_dim2 = ndpointer(dtype=np.int32, ndim=2, flags='C')

ipc_dll = CDLL("./build/Release/ipc.dll")
ipc_dll.init.argtypes = [nd_double_array_dim2, c_int, nd_int_array_dim2, c_int, nd_double_array_dim1, nd_double_array_dim1, nd_double_array_dim1, c_double, c_double]
ipc_dll.step.argtypes = []

n_seg = (4, 4)
side_len = 0.5
rho = 1000 # density of square
k = 1e4 # spring stiffness
dhat = 0.01
kappa = 1e5

x_length = (n_seg[0]+1)*(n_seg[1]+1)
e_length = n_seg[0]*(n_seg[1]+1) + (n_seg[0]+1)*n_seg[1] + n_seg[0]*n_seg[1]*2

x_data = np.zeros((x_length, 2), dtype=np.double)
e_data = np.zeros((e_length, 2), dtype=np.int32)
m_data = np.zeros(x_length, dtype=np.double)
l2_data = np.zeros(e_length, dtype=np.double)
k_data = np.zeros(e_length, dtype=np.double)

step = side_len / 4
for i in range(0, n_seg[0]+1):
    for j in range(0, n_seg[1]+1):
        x_data[i * (n_seg[1] + 1) + j] = [-side_len / 2.0 + i * step + 0.5, -side_len / 2.0 + j * step + 0.6]
e_i = 0
for i in range(0, n_seg[0]):
    for j in range(0, n_seg[1]+1):
        e_data[e_i] = [i * (n_seg[1] + 1) + j, (i + 1) * (n_seg[1] + 1) + j]
        e_i += 1
for i in range(0, n_seg[0]+1):
    for j in range(0, n_seg[1]):
        e_data[e_i] = [i * (n_seg[1] + 1) + j, i * (n_seg[1] + 1) + j + 1]
        e_i += 1
for i in range(0, n_seg[0]):
    for j in range(0, n_seg[1]):
        e_data[e_i] = [i * (n_seg[1] + 1) + j, (i + 1) * (n_seg[1] + 1) + j + 1]
        e_i += 1
        e_data[e_i] = [(i + 1) * (n_seg[1] + 1) + j, i * (n_seg[1] + 1) + j + 1]
        e_i += 1
for i in range(0, x_length):
    m_data[i] = rho * side_len * side_len / x_length
for i in range(0, e_length):
    k_data[i] = k
    edge = x_data[e_data[i][0]] - x_data[e_data[i][1]]
    l2_data[i] = edge[0] * edge[0] + edge[1] * edge[1]

ipc_dll.init(x_data, x_length, e_data, e_length, m_data, l2_data, k_data, dhat, kappa)

gui = ti.GUI("IPC")

while gui.running and not gui.get_event(gui.ESCAPE):
    ipc_dll.step()
    gui.clear(0xFFFFFF)
    gui.circles(x_data, radius=4, color=0x0000FF)
    for e in e_data:
        gui.line(x_data[e[0]], x_data[e[1]], radius=1, color=0x0000FF)
    gui.show()
