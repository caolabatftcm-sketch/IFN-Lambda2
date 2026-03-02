"""
Render IFN-Lambda alignment figures with variant sites annotated
(arrow + highlighted border).
Handles both IFNL2 (R154H) and IFNL3 (H132R, E175K).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math, os

# ── BLOSUM62 ──────────────────────────────────────────────────────────────────
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
    return BLOSUM62.get(key, BLOSUM62.get((b.upper(), a.upper()), -1))

def classify(h, m):
    if h == '-' or m == '-': return 'gap'
    if h.upper() == m.upper(): return 'identical'
    if blosum_score(h, m) > 0: return 'conservative'
    return 'noncons'

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


def build_pos_map(aln_h, aln_m):
    """
    Returns two dicts:
      h_to_col[human_1based_pos] = alignment_col_index (0-based)
      m_to_col[mouse_1based_pos] = alignment_col_index (0-based)
    And two reverse dicts:
      col_to_h[col] = human pos (or None)
      col_to_m[col] = mouse pos (or None)
    """
    h_to_col, m_to_col = {}, {}
    col_to_h, col_to_m = {}, {}
    h_pos = m_pos = 0
    for ci, (h, m) in enumerate(zip(aln_h, aln_m)):
        if h != '-':
            h_pos += 1
            h_to_col[h_pos] = ci
            col_to_h[ci] = h_pos
        else:
            col_to_h[ci] = None
        if m != '-':
            m_pos += 1
            m_to_col[m_pos] = ci
            col_to_m[ci] = m_pos
        else:
            col_to_m[ci] = None
    return h_to_col, m_to_col, col_to_h, col_to_m


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
VARIANT_BORDER = hex2rgb('C0392B')   # bold red border for variant cells
VARIANT_ARROW  = hex2rgb('C0392B')   # annotation arrow color


def render_alignment(aln_h, aln_m,
                     human_label, mouse_label,
                     title, subtitle, stats_txt,
                     variants,          # list of dicts: {h_pos, h_wt, h_mut, m_pos, m_aa, label}
                     out_path,
                     BLOCK=60, DPI=300):
    """
    variants: list of {
        'h_pos':  int,   human position (1-based)
        'h_wt':   str,   WT amino acid in human
        'h_mut':  str,   mutant amino acid
        'm_pos':  int or None,  corresponding mouse position
        'm_aa':   str,   mouse WT amino acid
        'label':  str,   e.g. 'R154H'
    }
    """
    aln_len = len(aln_h)

    # Layout
    FONT_SIZE  = 7
    LABEL_SIZE = 6.5
    CELL_W     = 0.135
    CELL_H     = 0.175
    CONS_H     = 0.11
    LABEL_W    = 1.30
    BLOCK_GAP  = 0.24
    ANNOT_H    = 0.38    # extra headroom per block that contains a variant annotation
    TOP_MARGIN = 0.58
    BOT_MARGIN = 0.20

    n_blocks = math.ceil(aln_len / BLOCK)

    # Which blocks contain variants?
    h_to_col, m_to_col, col_to_h, col_to_m = build_pos_map(aln_h, aln_m)
    variant_cols = set()
    for v in variants:
        if v['h_pos'] in h_to_col:
            variant_cols.add(h_to_col[v['h_pos']])

    def block_has_variant(b):
        s, e = b * BLOCK, min((b+1)*BLOCK, aln_len)
        return any(s <= c < e for c in variant_cols)

    # Compute per-block heights
    block_heights = []
    for b in range(n_blocks):
        bh = 2 * CELL_H + CONS_H
        if block_has_variant(b):
            bh += ANNOT_H
        block_heights.append(bh)

    fig_w = LABEL_W + BLOCK * CELL_W
    fig_h = (TOP_MARGIN + sum(block_heights) +
             (n_blocks - 1) * BLOCK_GAP + BOT_MARGIN)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title / subtitle / stats
    ax.text(fig_w / 2, fig_h - 0.08, title,
            ha='center', va='top', fontsize=9, fontweight='bold',
            fontfamily='Liberation Sans', color='#2E4057')
    ax.text(fig_w / 2, fig_h - 0.22, subtitle,
            ha='center', va='top', fontsize=6.5, color='#555555',
            fontfamily='Liberation Sans')
    ax.text(fig_w / 2, fig_h - 0.33, stats_txt,
            ha='center', va='top', fontsize=6.5, color='#333333',
            fontfamily='Liberation Sans')

    # Legend
    leg_y  = fig_h - 0.47
    leg_items = [
        ('Identical',        COLORS['identical']['bg'],    COLORS['identical']['fg']),
        ('Conservative',     COLORS['conservative']['bg'], COLORS['conservative']['fg']),
        ('Non-conservative', COLORS['noncons']['bg'],      COLORS['noncons']['fg']),
        ('Gap',              COLORS['gap']['bg'],           COLORS['gap']['fg']),
    ]
    leg_box_w, leg_box_h = 0.18, 0.10
    leg_x_start = 0.3
    leg_spacing = (fig_w - 0.6) / (len(leg_items) + 1)
    for i, (label, bg, fg) in enumerate(leg_items):
        lx = leg_x_start + i * leg_spacing
        rect = mpatches.FancyBboxPatch((lx, leg_y - leg_box_h/2), leg_box_w, leg_box_h,
                                       boxstyle="round,pad=0.01", linewidth=0.3,
                                       edgecolor='#aaaaaa', facecolor=bg)
        ax.add_patch(rect)
        ax.text(lx + leg_box_w + 0.04, leg_y, label,
                va='center', fontsize=6, fontfamily='Liberation Sans', color='#333333')
    # Variant marker in legend
    lx = leg_x_start + len(leg_items) * leg_spacing
    rect_v = mpatches.FancyBboxPatch((lx, leg_y - leg_box_h/2), leg_box_w, leg_box_h,
                                      boxstyle="round,pad=0.01", linewidth=1.5,
                                      edgecolor=VARIANT_BORDER, facecolor=hex2rgb('fdecea'))
    ax.add_patch(rect_v)
    ax.text(lx + leg_box_w + 0.04, leg_y, 'Variant site',
            va='center', fontsize=6, fontfamily='Liberation Sans', color='#333333')

    ax.text(fig_w - 1.05, leg_y, '| = identical   + = conservative',
            va='center', fontsize=5.8, fontfamily='Liberation Sans', color='#555555')

    # Build variant col → info lookup
    vcol_info = {}
    for v in variants:
        if v['h_pos'] in h_to_col:
            vcol_info[h_to_col[v['h_pos']]] = v

    # Draw blocks
    y_cursor = fig_h - TOP_MARGIN
    for b in range(n_blocks):
        s = b * BLOCK
        e = min(s + BLOCK, aln_len)
        bh_seq = aln_h[s:e]
        bm_seq = aln_m[s:e]
        has_var = block_has_variant(b)
        annot_h = ANNOT_H if has_var else 0

        if b > 0:
            ax.plot([0.05, fig_w - 0.05],
                    [y_cursor + 0.015, y_cursor + 0.015],
                    color='#dddddd', linewidth=0.4)

        # Row tops (accounting for annotation headroom at top of block)
        y_human = y_cursor - annot_h
        y_cons  = y_human - CELL_H
        y_mouse = y_cons  - CONS_H

        h_start = sum(1 for c in aln_h[:s] if c != '-') + 1
        h_end   = sum(1 for c in aln_h[:e] if c != '-')
        m_start = sum(1 for c in aln_m[:s] if c != '-') + 1
        m_end   = sum(1 for c in aln_m[:e] if c != '-')

        # Labels
        ax.text(LABEL_W - 0.06, y_human - CELL_H/2,
                f'{human_label} {h_start}–{h_end}', ha='right', va='center',
                fontsize=LABEL_SIZE, fontfamily='DejaVu Sans Mono',
                color='#2E4057', fontweight='bold')
        ax.text(LABEL_W - 0.06, y_mouse - CELL_H/2,
                f'{mouse_label} {m_start}–{m_end}', ha='right', va='center',
                fontsize=LABEL_SIZE, fontfamily='DejaVu Sans Mono',
                color='#2E4057', fontweight='bold')

        # Per-column residues
        for ci, (h, m) in enumerate(zip(bh_seq, bm_seq)):
            global_col = s + ci
            cls = classify(h, m)
            bg  = COLORS[cls]['bg']
            fg  = COLORS[cls]['fg']
            x   = LABEL_W + ci * CELL_W
            is_variant = global_col in vcol_info

            # Human cell
            rect_h = mpatches.Rectangle(
                (x, y_human - CELL_H), CELL_W, CELL_H,
                facecolor=hex2rgb('fdecea') if is_variant else bg,
                edgecolor=VARIANT_BORDER if is_variant else 'white',
                linewidth=1.4 if is_variant else 0.25, zorder=3 if is_variant else 1)
            ax.add_patch(rect_h)
            ax.text(x + CELL_W/2, y_human - CELL_H/2, h,
                    ha='center', va='center', fontsize=FONT_SIZE,
                    fontfamily='DejaVu Sans Mono',
                    color=VARIANT_BORDER if is_variant else fg,
                    fontweight='bold', zorder=4)

            # Mouse cell
            rect_m = mpatches.Rectangle(
                (x, y_mouse - CELL_H), CELL_W, CELL_H,
                facecolor=hex2rgb('fdecea') if is_variant else bg,
                edgecolor=VARIANT_BORDER if is_variant else 'white',
                linewidth=1.4 if is_variant else 0.25, zorder=3 if is_variant else 1)
            ax.add_patch(rect_m)
            ax.text(x + CELL_W/2, y_mouse - CELL_H/2, m,
                    ha='center', va='center', fontsize=FONT_SIZE,
                    fontfamily='DejaVu Sans Mono',
                    color=VARIANT_BORDER if is_variant else fg,
                    fontweight='bold', zorder=4)

            # Consensus
            if cls == 'identical':
                ax.text(x + CELL_W/2, y_cons - CONS_H/2, '|',
                        ha='center', va='center', fontsize=FONT_SIZE - 0.5,
                        fontfamily='DejaVu Sans Mono',
                        color=CONS_COLORS['identical'], fontweight='bold')
            elif cls == 'conservative':
                ax.text(x + CELL_W/2, y_cons - CONS_H/2, '+',
                        ha='center', va='center', fontsize=FONT_SIZE - 0.5,
                        fontfamily='DejaVu Sans Mono',
                        color=CONS_COLORS['conservative'], fontweight='bold')

            # Variant annotation: single label line + thin tick above human cell
            if is_variant and has_var:
                v = vcol_info[global_col]
                cell_cx  = x + CELL_W / 2
                cell_top = y_human              # top edge of human cell
                tick_bot = cell_top + 0.01
                tick_top = cell_top + ANNOT_H * 0.45   # bottom of text area
                label_y  = cell_top + ANNOT_H * 0.92   # center of label

                # Single label: "H132R  →  Mm H124"
                m_note = f'{v["m_aa"]}{v["m_pos"]}' if v['m_pos'] else 'gap'
                full_label = f'{v["h_wt"]}{v["h_pos"]}{v["h_mut"]}  \u2192  Mm {m_note}'
                ax.text(cell_cx, label_y, full_label,
                        ha='center', va='center',
                        fontsize=6.2, fontweight='bold',
                        fontfamily='Liberation Sans',
                        color=VARIANT_BORDER, zorder=10,
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                  edgecolor=VARIANT_BORDER, linewidth=0.8, alpha=0.92))

                # Thin vertical tick from label to cell top
                ax.plot([cell_cx, cell_cx], [tick_bot, tick_top],
                        color=VARIANT_BORDER, lw=1.0, zorder=9)
                # Small downward arrowhead at bottom of tick
                ax.annotate('',
                    xy=(cell_cx, tick_bot),
                    xytext=(cell_cx, tick_top),
                    arrowprops=dict(arrowstyle='->', color=VARIANT_BORDER,
                                    lw=0.9, mutation_scale=7),
                    zorder=10)

        y_cursor -= (block_heights[b] + BLOCK_GAP)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# 1. IFN-λ2: Human vs Mouse  —  R154H
# ══════════════════════════════════════════════════════════════════════════════
HUMAN_L2 = "MKLDMTGDCTPVLVLMAAVLTVTGAVPVARLHGALPDARGCHIAQFKSLSPQELQAFKRAKDALEESLLLKDCRCHSRLFPRTWDLRQLQVRERPMALEAELALTLKVLEATADTDPALVDVLDQPLHTLHHILSQFRACIQPQPTAGPRTRGRLHHWLYRLQEAPKKESPGCLEASVTFNLFRLLTRDLNCVASGDLCV"
MOUSE_L2  = "MLLLLLPLLLAAVLTRTQADPVPRATRLPVEAKDCHIAQFKSLSPKELQAFKKAKDAIEKRLLEKDLRCSSHLFPRAWDLKQLQVQERPKALQAEVALTLKVWENMTDSALATILGQPLHTLSHIHSQLQTCTQLQATAEPRSPSRRLSRWLHRLQEAQSKETPGCLEASVTSNLFRLLTRDLKCVANGDQCV"

print("Aligning IFNL2...")
aln_h2, aln_m2 = needleman_wunsch(HUMAN_L2, MOUSE_L2)
h2_to_col, m2_to_col, col2_to_h, col2_to_m = build_pos_map(aln_h2, aln_m2)

# Find mouse position for human R154
col_154 = h2_to_col.get(154)
m_pos_154 = col2_to_m.get(col_154) if col_154 is not None else None
m_aa_154  = aln_m2[col_154] if col_154 is not None else '?'
print(f"  Human R154 → alignment col {col_154} → Mouse pos {m_pos_154} ({m_aa_154})")

paired2 = [(aln_h2[i], aln_m2[i]) for i in range(len(aln_h2))
           if aln_h2[i] != '-' and aln_m2[i] != '-']
ident2  = sum(1 for h, m in paired2 if h.upper() == m.upper())
cons2   = sum(1 for h, m in paired2 if h != m and blosum_score(h,m) > 0)
pct_id2  = 100 * ident2 / len(paired2)
pct_sim2 = 100 * (ident2 + cons2) / len(paired2)

variants_l2 = [{
    'h_pos': 154, 'h_wt': 'R', 'h_mut': 'H',
    'm_pos': m_pos_154, 'm_aa': m_aa_154,
    'label': 'R154H',
}]

render_alignment(
    aln_h2, aln_m2,
    human_label='Hs', mouse_label='Mm',
    title='IFN-Lambda 2: Human (Q8IZJ0) vs Mouse (Q4VK74) — Variant R154H annotated',
    subtitle='Needleman–Wunsch global alignment  ·  BLOSUM62  ·  Gap open –10  ·  Gap extend –0.5',
    stats_txt=(f'Identity: {ident2}/{len(paired2)} ({pct_id2:.1f}%)   '
               f'Similarity: {ident2+cons2}/{len(paired2)} ({pct_sim2:.1f}%)   '
               f'Alignment length: {len(aln_h2)}'),
    variants=variants_l2,
    out_path='/sessions/affectionate-zealous-hypatia/mnt/outputs/IFNL2_human_mouse_R154H_annotated.png',
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. IFN-λ3: Human vs Mouse  —  H132R, E175K
# ══════════════════════════════════════════════════════════════════════════════
# Human IFNL3 WT isoform 1 (from document)
HUMAN_L3 = "MKLDMTGDCMPVLVLMAAVLTVTGAVPVARLRGALPDARGCHIAQFKSLSPQELQAFKRAKDALEESLLLKDCKCRSRLFPRTWDLRQLQVRERPVALEAELALTLKVLEATADTDPALGDVLDQPLHTLHHILSQLRACIQPQPTAGPRTRGRLHHWLHRLQEAPKKESPGCLEASVTFNLFRLLTRDLNCVASGDLCV"
MOUSE_L3  = "MLLLLLPLLLAAVLTRTQADPVPRATRLPVEAKDCHIAQFKSLSPKELQAFKKAKGAIEKRLLEKDMRCSSHLISRAWDLKQLQVQERPKALQAEVALTLKVWENINDSALTTILGQPLHTLSHIHSQLQTCTQLQATAEPKPPSRRLSRWLHRLQEAQSKETPGCLEDSVTSNLFQLLLRDLKCVASGDQCV"

print("Aligning IFNL3...")
aln_h3, aln_m3 = needleman_wunsch(HUMAN_L3, MOUSE_L3)
h3_to_col, m3_to_col, col3_to_h, col3_to_m = build_pos_map(aln_h3, aln_m3)

for hpos, hlabel in [(132, 'H132'), (175, 'E175')]:
    col = h3_to_col.get(hpos)
    mpos = col3_to_m.get(col) if col is not None else None
    maa  = aln_m3[col] if col is not None else '?'
    print(f"  Human {hlabel} → alignment col {col} → Mouse pos {mpos} ({maa})")

col_132 = h3_to_col.get(132); m_pos_132 = col3_to_m.get(col_132); m_aa_132 = aln_m3[col_132] if col_132 else '?'
col_175 = h3_to_col.get(175); m_pos_175 = col3_to_m.get(col_175); m_aa_175 = aln_m3[col_175] if col_175 else '?'

paired3 = [(aln_h3[i], aln_m3[i]) for i in range(len(aln_h3))
           if aln_h3[i] != '-' and aln_m3[i] != '-']
ident3  = sum(1 for h, m in paired3 if h.upper() == m.upper())
cons3   = sum(1 for h, m in paired3 if h != m and blosum_score(h,m) > 0)
pct_id3  = 100 * ident3 / len(paired3)
pct_sim3 = 100 * (ident3 + cons3) / len(paired3)

variants_l3 = [
    {'h_pos': 132, 'h_wt': 'H', 'h_mut': 'R',
     'm_pos': m_pos_132, 'm_aa': m_aa_132, 'label': 'H132R'},
    {'h_pos': 175, 'h_wt': 'E', 'h_mut': 'K',
     'm_pos': m_pos_175, 'm_aa': m_aa_175, 'label': 'E175K'},
]

render_alignment(
    aln_h3, aln_m3,
    human_label='Hs', mouse_label='Mm',
    title='IFN-Lambda 3: Human (isoform 1) vs Mouse (Q8CGK6) — Variants H132R, E175K annotated',
    subtitle='Needleman–Wunsch global alignment  ·  BLOSUM62  ·  Gap open –10  ·  Gap extend –0.5',
    stats_txt=(f'Identity: {ident3}/{len(paired3)} ({pct_id3:.1f}%)   '
               f'Similarity: {ident3+cons3}/{len(paired3)} ({pct_sim3:.1f}%)   '
               f'Alignment length: {len(aln_h3)}'),
    variants=variants_l3,
    out_path='/sessions/affectionate-zealous-hypatia/mnt/outputs/IFNL3_human_mouse_H132R_E175K_annotated.png',
)
print("All done.")
