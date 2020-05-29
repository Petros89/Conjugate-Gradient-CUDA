import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

labels = ['25', '50', '100', '150', '200']
mat_free = [2, 5, 10, 18, 30]
mat_glob = [5, 12, 24, 41, 63]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, mat_free, width, label='Matrix-Free')
rects2 = ax.bar(x + width/2, mat_glob, width, label='Global-Matrix')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Wall Time [seconds]')
ax.set_xlabel('Mesh Size [$element^3$]')
ax.set_title('Matrix Free VS Global Matrix Timings [CPU]')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

plt.show()

