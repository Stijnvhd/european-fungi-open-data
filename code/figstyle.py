import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
VIOLET = "#6250d6"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
SEQ = ["#d6e6f8", "#a9cbf0", "#6fa7e4", "#2a78d6", "#1a56a1", "#0f3768"]

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300,
    "font.family": "DejaVu Sans", "font.size": 9.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.size": 0, "ytick.major.size": 0,
    "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.titlelocation": "left", "axes.titlepad": 10,
    "legend.frameon": False, "legend.fontsize": 8.5,
    "lines.linewidth": 2.0,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def save(fig, name):
    fig.savefig(f"figures/{name}.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def source(ax, text):
    ax.annotate(text, xy=(0, -0.16), xycoords="axes fraction", fontsize=7, color=MUTED, ha="left", va="top")
