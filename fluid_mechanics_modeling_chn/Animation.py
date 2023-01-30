from __future__ import division
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

from main_program import *


# 定义读取模拟中数据的函数
# read the data from simulation
def read_datafile(iteration):
    filename = "PUV{0}.txt".format(iteration)
    filepath = os.path.join(dir_path, filename)
    arr = np.loadtxt(filepath, delimiter="\t")
    rows, cols = arr.shape
    p_p = np.zeros((rowpts, colpts))
    u_p = np.zeros((rowpts, colpts))
    v_p = np.zeros((rowpts, colpts))
    p_arr = arr[:, 0]
    u_arr = arr[:, 1]
    v_arr = arr[:, 2]

    p_p = p_arr.reshape((rowpts, colpts))
    u_p = u_arr.reshape((rowpts, colpts))
    v_p = v_arr.reshape((rowpts, colpts))

    return p_p, u_p, v_p


# 结果文件
# the result file
cwdir = os.getcwd()
dir_path = os.path.join(cwdir, "Result")
os.chdir(dir_path)

# 遍历directory中的文件
# go through all files in the directory
filenames = []
iterations = []
for root, dirs, files in os.walk(dir_path):
    for datafile in files:
        if "PUV" in datafile:
            filenames.append(datafile)
            no_ext_file = datafile.replace(".txt", "").strip()
            iter_no = int(no_ext_file.split("V")[-1])
            iterations.append(iter_no)

# 辨别最后的迭代与迭代数
# identify the iteration number and the last iteration
initial_iter = np.amin(iterations)
final_iter = np.amax(iterations)
inter = (final_iter - initial_iter) / (len(iterations) - 1)
number_of_frames = len(iterations)  # int(final_iter/inter)+1
sorted_iterations = np.sort(iterations)

# 创建数组和网格
# set up array and grids
x = np.linspace(0, length, colpts)
y = np.linspace(0, breadth, rowpts)
[X, Y] = np.meshgrid(x, y)

# 定义streamplot的编号
# define the index of streamplot
index_cut_x = int(colpts / 10)
index_cut_y = int(rowpts / 10)

# 产生空白图片
# generate a blank picture
fig = plt.figure(figsize=(16, 8))
ax = plt.axes(xlim=(0, length), ylim=(0, breadth))

# 初始等高线图
# initialize the contour plot
p_p, u_p, v_p = read_datafile(0)
ax.set_xlim([0, length])
ax.set_ylim([0, breadth])
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_title("Frame No: 0")
cont = ax.contourf(X, Y, p_p)
stream = ax.streamplot(X[::index_cut_y, ::index_cut_x], Y[::index_cut_y, ::index_cut_x],
                       u_p[::index_cut_y, ::index_cut_x], v_p[::index_cut_y, ::index_cut_x], color="k")
fig.colorbar(cont)
fig.tight_layout()


def animate(i):
    sys.stdout.write("\rframe remain: {0:03d}".format(len(sorted_iterations) - i))
    sys.stdout.flush()
    iteration = sorted_iterations[i]
    p_p, u_p, v_p = read_datafile(iteration)
    ax.clear()
    ax.set_xlim([0, length])
    ax.set_ylim([0, breadth])
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel("$y$", fontsize=12)
    ax.set_title("Frame No: {0}".format(i))
    cont = ax.contourf(X, Y, p_p)
    stream = ax.streamplot(X[::index_cut_y, ::index_cut_x], Y[::index_cut_y, ::index_cut_x],
                           u_p[::index_cut_y, ::index_cut_x], v_p[::index_cut_y, ::index_cut_x], color="k")
    return cont, stream



