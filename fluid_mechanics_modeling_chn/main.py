from Animation import *
from head_file import *
from main_program import *

if __name__ == '__main__':
    print("########        animation        ########")
    print("#########################################")
    anim = animation.FuncAnimation(fig, animate, frames=number_of_frames, interval=50, blit=False)

    # save the simulation gif to movie_path with PillowWriter
    movie_path = os.path.join(dir_path, "animation.gif")
    anim.save(r"{0}".format(movie_path), writer=PillowWriter())
    print("\n动画在Result中保存为animation.gif")
    print("\nsaved as animation.gif in Result")
