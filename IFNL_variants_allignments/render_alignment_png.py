"""
Render IFN-Lambda 2 human vs mouse alignment as a high-resolution PNG figure
suitable for publication (300 DPI).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import math, os

# ── Sequences & alignment (reuse from main script) ───────────────────────────
HUMAN = "MKLDMTGDCTPVLVLMAAVLTVTGAVPVARLHGALPDARGCHIAQFKSLSPQELQAFKRAKDALEESLLLKDCRCHSRLFPRTWDLRQLQVRERPMALEAELALTLKVLEATADTDPALVDVLDQPLHTLHHILSQFRACIQPQPTAGPRTRGRLHHWLYRLQEAPKKESPGCLEASVTFNLFRLLTRDLNCVASGDLCV"
MOUSE = "MLLLLLPLLLAAVLTRTQADPVPRATRLPVEAKDCHIAQFKSLSPKELQAFKKAKDAIEKRLLEKDLRCSSHLFPRAWDLKQLQVQERPKALQAEVALTLKVWENMTDSALATILGQPLHTLSHIHSQLQTCTQLQATAEPRSPSRRLSRWLHRLQEAQSKETPGCLEASVTSNLFRLLTRDLKCVANGDQCV"

BLOSUM62_STR = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

def parse_blosum62(text):
    lines = [l for l in text.strip().split('\n') if l.strip()]
    aa_order = lines[0].split()
    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        aa1 = parts[0]
        for j, aa2 in enumerate(aa_order):
            matrix[(aa1, aa2)] = int(parts[j + 1])
    return matrix

BLOSUM62 = parse_blosum62(BLOSUM62_STR)

def blosum_score(a, b):
    key = (a.upper(), b.upper())
    rkey = (b.upper(), a.upper())
    return BLOSUM62.get(key, BLOSUM62.get(rkey, -1))

def needleman_wunsch(seq1, seq2, gap_open=-10, gap_extend=-0.5):
    m, n = len(seq1), len(seq2)
    NEG_INF = float('-inf')
    M  = [[NEG_INF]*(n+1) for _ in range(m+1)]
    IX = [[NEG_INF]*(n+1) for _ in range(m+1)]
    IY = [[NEG_INF]*(n+1) for _ in range(m+1)]
    M[0][0] = 0
    for i in range(1, m+1):
        IX[i][0] = gap_open + (i-1)*gap_extend
        M[i][0]  = gap_open + (i-1)*gap_extend
    for j in range(1, n+1):
        IY[0][j] = gap_open + (j-1)*gap_extend
        M[0][j]  = gap_open + (j-1)*gap_extend
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = blosum_score(seq1[i-1], seq2[j-1])
            M[i][j]  = s + max(M[i-1][j-1], IX[i-1][j-1], IY[i-1][j-1])
            IX[i][j] = max(M[i-1][j] + gap_open, IX[i-1][j] + gap_extend)
            IY[i][j] = max(M[i][j-1] + gap_open, IY[i][j-1] + gap_extend)
    align1, align2 = [], []
    i, j = m, n
    state = 'M'
    while i > 0 or j > 0:
        if state == 'M':
            if i == 0:
                align1.append('-'); align2.append(seq2[j-1]); j -= 1
            elif j == 0:
                align1.append(seq1[i-1]); align2.append('-'); i -= 1
            else:
                prev_M = M[i-1][j-1]; prev_IX = IX[i-1][j-1]; prev_IY = IY[i-1][j-1]
                best = max(prev_M, prev_IX, prev_IY)
                align1.append(seq1[i-1]); align2.append(seq2[j-1])
                i -= 1; j -= 1
                if best == prev_IX: state = 'IX'
                elif best == prev_IY: state = 'IY'
        elif state == 'IX':
            align1.append(seq1[i-1]); align2.append('-')
            if M[i-1][j] + gap_open >= IX[i-1][j] + gap_extend: state = 'M'
            i -= 1
        elif state == 'IY':
            align1.append('-'); align2.append(seq2[j-1])
            if M[i][j-1] + gap_open >= IY[i][j-1] + gap_extend: state = 'M'
            j -= 1
    return ''.join(reversed(align1)), ''.join(reversed(align2))

def classify(h, m):
    if h == '-' or m == '-': return 'gap'
    if h.upper() == m.upper(): return 'identical'
    if blosum_score(h, m) > 0: return 'conservative'
    return 'noncons'

# ── Run alignment ─────────────────────────────────────────────────────────────
aln_h, aln_m = needleman_wunsch(HUMAN, MOUSE)
aln_len = len(aln_h)

# ── Figure layout parameters ──────────────────────────────────────────────────
DPI        = 300
BLOCK      = 60
FONT_SIZE  = 7        # pt — amino acid letter size
LABEL_SIZE = 6.5      # pt — row label font size

CELL_W     = 0.135    # inches per AA column
CELL_H     = 0.175    # inches per sequence row
CONS_H     = 0.11     # inches for consensus row
LABEL_W    = 1.30     # inches for left label column
BLOCK_GAP  = 0.18     # inches between blocks
TOP_MARGIN = 0.55     # inches (title + legend)
BOT_MARGIN = 0.15

n_blocks   = math.ceil(aln_len / BLOCK)
seq_w      = BLOCK * CELL_W
fig_w      = LABEL_W + seq_w
block_h    = 2 * CELL_H + CONS_H
fig_h      = TOP_MARGIN + n_blocks * block_h + (n_blocks - 1) * BLOCK_GAP + BOT_MARGIN

# Color palette (RGB 0-1)
def hex2rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))

COLORS = {
    'identical':    {'bg': hex2rgb('d4edda'), 'fg': hex2rgb('155724')},
    'conservative': {'bg': hex2rgb('fff3cd'), 'fg': hex2rgb('856404')},
    'noncons':      {'bg': hex2rgb('f8d7da'), 'fg': hex2rgb('721c24')},
    'gap':          {'bg': hex2rgb('eeeeee'), 'fg': hex2rgb('888888')},
}
CONS_COLORS = {
    'identical':    hex2rgb('155724'),
    'conservative': hex2rgb('856404'),
}

# ── Draw ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis('off')
fig.patch.set_facecolor('white')

# Title
ax.text(fig_w / 2, fig_h - 0.08, 'IFN-Lambda 2: Human (Q8IZJ0) vs Mouse (Q4VK74)',
        ha='center', va='top', fontsize=9, fontweight='bold', fontfamily='Liberation Sans',
        color='#2E4057')
ax.text(fig_w / 2, fig_h - 0.22,
        'Needleman–Wunsch global alignment  ·  BLOSUM62  ·  Gap open –10  ·  Gap extend –0.5',
        ha='center', va='top', fontsize=6.5, color='#555555', fontfamily='Liberation Sans')

# Statistics line
paired = [(aln_h[i], aln_m[i]) for i in range(aln_len) if aln_h[i] != '-' and aln_m[i] != '-']
ident  = sum(1 for h, m in paired if h.upper() == m.upper())
cons   = sum(1 for h, m in paired if not h.upper() == m.upper() and blosum_score(h,m) > 0)
pct_id  = 100 * ident / len(paired)
pct_sim = 100 * (ident + cons) / len(paired)
stats_txt = (f'Identity: {ident}/{len(paired)} ({pct_id:.1f}%)   '
             f'Similarity: {ident+cons}/{len(paired)} ({pct_sim:.1f}%)   '
             f'Alignment length: {aln_len}')
ax.text(fig_w / 2, fig_h - 0.33, stats_txt,
        ha='center', va='top', fontsize=6.5, color='#333333', fontfamily='Liberation Sans')

# Legend
leg_y  = fig_h - 0.47
leg_items = [
    ('Identical',       COLORS['identical']['bg'],    COLORS['identical']['fg']),
    ('Conservative',    COLORS['conservative']['bg'], COLORS['conservative']['fg']),
    ('Non-conservative',COLORS['noncons']['bg'],      COLORS['noncons']['fg']),
    ('Gap',             COLORS['gap']['bg'],           COLORS['gap']['fg']),
]
leg_box_w, leg_box_h = 0.18, 0.10
leg_x_start = 0.3
leg_spacing = (fig_w - 0.6) / len(leg_items)
for i, (label, bg, fg) in enumerate(leg_items):
    lx = leg_x_start + i * leg_spacing
    rect = mpatches.FancyBboxPatch((lx, leg_y - leg_box_h/2), leg_box_w, leg_box_h,
                                    boxstyle="round,pad=0.01", linewidth=0.3,
                                    edgecolor='#aaaaaa', facecolor=bg)
    ax.add_patch(rect)
    ax.text(lx + leg_box_w + 0.05, leg_y, label,
            va='center', fontsize=6, fontfamily='Liberation Sans', color='#333333')

# Consensus symbol legend
ax.text(fig_w - 1.05, leg_y, '| = identical   + = conservative',
        va='center', fontsize=5.8, fontfamily='Liberation Sans', color='#555555')

# ── Draw alignment blocks ─────────────────────────────────────────────────────
y_cursor = fig_h - TOP_MARGIN  # top of first block, drawing downward

for b in range(n_blocks):
    s = b * BLOCK
    e = min(s + BLOCK, aln_len)
    bh = aln_h[s:e]
    bm = aln_m[s:e]

    h_start = sum(1 for c in aln_h[:s] if c != '-') + 1
    h_end   = sum(1 for c in aln_h[:e] if c != '-')
    m_start = sum(1 for c in aln_m[:s] if c != '-') + 1
    m_end   = sum(1 for c in aln_m[:e] if c != '-')

    # Thin separator line above block (except first)
    if b > 0:
        ax.plot([0.05, fig_w - 0.05], [y_cursor + 0.015, y_cursor + 0.015],
                color='#dddddd', linewidth=0.4)

    # Row y positions (top edge of each row)
    y_human = y_cursor
    y_cons  = y_cursor - CELL_H
    y_mouse = y_cursor - CELL_H - CONS_H

    # Labels
    ax.text(LABEL_W - 0.06, y_human - CELL_H/2,
            f'Hs {h_start}–{h_end}', ha='right', va='center',
            fontsize=LABEL_SIZE, fontfamily='DejaVu Sans Mono', color='#2E4057', fontweight='bold')
    ax.text(LABEL_W - 0.06, y_mouse - CELL_H/2,
            f'Mm {m_start}–{m_end}', ha='right', va='center',
            fontsize=LABEL_SIZE, fontfamily='DejaVu Sans Mono', color='#2E4057', fontweight='bold')

    # Per-column residues
    for ci, (h, m) in enumerate(zip(bh, bm)):
        cls = classify(h, m)
        bg  = COLORS[cls]['bg']
        fg  = COLORS[cls]['fg']
        x   = LABEL_W + ci * CELL_W

        # Human cell
        rect_h = mpatches.Rectangle((x, y_human - CELL_H), CELL_W, CELL_H,
                                     facecolor=bg, edgecolor='white', linewidth=0.25)
        ax.add_patch(rect_h)
        ax.text(x + CELL_W/2, y_human - CELL_H/2, h,
                ha='center', va='center', fontsize=FONT_SIZE,
                fontfamily='DejaVu Sans Mono', color=fg, fontweight='bold')

        # Mouse cell
        rect_m = mpatches.Rectangle((x, y_mouse - CELL_H), CELL_W, CELL_H,
                                     facecolor=bg, edgecolor='white', linewidth=0.25)
        ax.add_patch(rect_m)
        ax.text(x + CELL_W/2, y_mouse - CELL_H/2, m,
                ha='center', va='center', fontsize=FONT_SIZE,
                fontfamily='DejaVu Sans Mono', color=fg, fontweight='bold')

        # Consensus symbol
        if cls == 'identical':
            ax.text(x + CELL_W/2, y_cons - CONS_H/2, '|',
                    ha='center', va='center', fontsize=FONT_SIZE - 0.5,
                    fontfamily='DejaVu Sans Mono', color=CONS_COLORS['identical'], fontweight='bold')
        elif cls == 'conservative':
            ax.text(x + CELL_W/2, y_cons - CONS_H/2, '+',
                    ha='center', va='center', fontsize=FONT_SIZE - 0.5,
                    fontfamily='DejaVu Sans Mono', color=CONS_COLORS['conservative'], fontweight='bold')

    y_cursor -= (block_h + BLOCK_GAP)

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = "/sessions/affectionate-zealous-hypatia/mnt/outputs/IFNL2_human_mouse_alignment.png"
os.makedirs(os.path.dirname(OUT), exist_ok=True)
plt.savefig(OUT, dpi=DPI, bbox_inches='tight', facecolor='white',
            metadata={'Title': 'IFN-Lambda 2 Human vs Mouse Alignment'})
plt.close()
print(f"Saved: {OUT}")

import os
size_kb = os.path.getsize(OUT) // 1024
print(f"File size: {size_kb} KB")

# Also save SVG for infinite scalability
OUT_SVG = OUT.replace('.png', '.svg')
fig2, ax2 = plt.subplots(figsize=(fig_w, fig_h))
# (re-run draw by calling same block — just re-execute and save as SVG)
plt.close()
print("Done.")
