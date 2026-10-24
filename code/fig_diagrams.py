import sys
sys.path.insert(0, "scripts")
from figstyle import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def box(ax, x, y, w, h, title, lines, fc="#eef3fa", ec=BLUE, fs=7.2, left=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h - 0.06, title, ha="center", va="top", fontsize=fs, color=INK, weight="bold")
    body = "\n".join(("• " + l) if left else l for l in lines)
    ax.text(x + (0.025 if left else w / 2), y + h - 0.2, body, ha="left" if left else "center", va="top",
            fontsize=fs - 0.6, color=INK2, linespacing=1.5)


def arrow(ax, x0, y0, x1, y1, col=MUTED):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, color=col, lw=1.2))


def blank(figsize):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    return fig, ax


def flow(name, steps, note, fc, ec):
    fig, ax = blank((6.8, 2.2))
    w, h, gap = 0.17, 0.5, 0.035
    for i, (t, ls) in enumerate(steps):
        x = 0.01 + i * (w + gap)
        box(ax, x, 0.4, w, h, t, ls, fc=fc, ec=ec)
        if i < len(steps) - 1:
            arrow(ax, x + w + 0.004, 0.65, x + w + gap - 0.004, 0.65)
    ax.text(0.01, 0.03, note, fontsize=6.6, color=MUTED, va="bottom")
    save(fig, name)


flow("record_journey",
     [("In the wood", ["a recorder finds", "and names", "a fungus"]), ("The scheme", ["a society or", "platform checks", "and stores it"]),
      ("GBIF", ["common format,", "names matched", "to the backbone"]), ("The query", ["counts by", "species, place", "and ten days"]),
      ("The analysis", ["effort, scheme", "by scheme, tests", "of alternatives"])],
     "Where bias enters: who goes out, when and where (effort); how names are chosen (identification);\nwhich schemes publish openly (coverage); how names are matched (taxonomy).", "#eef3fa", BLUE)

flow("metabarcoding_steps",
     [("Soil", ["a quarter gram;", "DNA of all life", "extracted"]), ("Copy", ["primers flank", "ITS; PCR makes", "many copies"]),
      ("Sequence", ["tens of thousands", "of tagged reads", "per sample"]), ("Clean", ["errors removed;", "similar reads", "clustered"]),
      ("Match", ["against UNITE:", "named, code only", "or unplaced"])],
     "Each step shapes the result: primer bias, copy-number differences between species, patchy sampling,\nand above all the reference library.", "#eef8f3", AQUA)

fig, ax = blank((6.8, 2.5))
box(ax, 0.01, 0.2, 0.47, 0.74, "Living with trees (ectomycorrhizal)",
    ["sheathes the root tips of a living tree", "receives sugar made in the leaves", "passes nitrogen and phosphorus back",
     "fruits on current and stored host sugar", "sensitive to extra nitrogen from the air", "chanterelle, cep, fly agaric, webcaps"], fs=7.6, left=True)
box(ax, 0.52, 0.2, 0.47, 0.74, "Rotting wood and leaves (saprotrophic)",
    ["digests dead leaves, needles and wood", "its energy is already in the substrate", "no partner to pay or be paid by",
     "fruits when temperature and moisture allow", "no tie to the sugar of living trees", "sulphur tuft, clouded agaric, velvet shank"],
    fc="#fdf0ea", ec=ORANGE, fs=7.6, left=True)
save(fig, "two_guilds")
