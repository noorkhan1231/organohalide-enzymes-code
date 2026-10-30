"""Shared figure style: Liberation Sans (Arial-metric), 6.5-7 pt, host colours used in every manuscript figure."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5, 'axes.labelweight': 'bold',
                     'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 4.0, 'ytick.major.size': 4.0, 'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0, 'xtick.minor.width': 0.45, 'ytick.minor.width': 0.45, 'savefig.dpi': 600,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic', 'mathtext.bf': 'Liberation Sans:bold'})
import matplotlib.legend as _ml
_orig_legend_init = _ml.Legend.__init__
def _bold_legend_init(self, *a, **k):
    _orig_legend_init(self, *a, **k)
    for t in self.get_texts(): t.set_fontweight('bold')
    if self.get_title() is not None: self.get_title().set_fontweight('bold')
_ml.Legend.__init__ = _bold_legend_init
HOSTS = ['Chicken', 'Coral', 'Rhizosphere', 'Sponge']
COL = {'Chicken': '#E69F00', 'Coral': '#56B4E9', 'Rhizosphere': '#009E73', 'Sponge': '#CC79A7'}
FAM = ['RDase', 'HyDase', 'GSTDase', 'OxDase', 'DehyDase', 'MetDase', 'FlaHase', 'VanHase', 'NHFeHase', 'SAMHase']
def label(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s.upper(), transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom', ha='left')
def save(fig, out):
    fig.savefig(out + '.pdf'); fig.savefig(out + '.jpg', dpi=1200, pil_kwargs={'quality': 95}); fig.savefig(out + '_preview.png', dpi=200)
def pfmt(p):
    if p < 0.001: return 'P < 0.001'
    return f'P = {p:.3f}'

MM = 1 / 25.4          # inches per mm
FULL_W = 174           # Cell Press full-page width (mm)

def assemble(panels, out, width_mm=FULL_W, letters=True, dpi=1200):
    """Place panel PDFs (vector) on one page. panels: list of (letter, pdf_path, x_mm, y_mm, w_mm, h_mm).
    Each panel is scaled to fit its box (aspect kept, top-left aligned). Writes out.pdf, out.jpg (dpi), out_preview.png."""
    import pymupdf
    pt = 72 / 25.4
    H = max(y + h for _, _, x, y, w, h in panels) + 1.0
    doc = pymupdf.open(); page = doc.new_page(width=width_mm * pt, height=H * pt)
    for let, f, x, y, w, h in panels:
        src = pymupdf.open(f); r = src[0].rect; s = min(w * pt / r.width, h * pt / r.height)
        rect = pymupdf.Rect(x * pt, y * pt, x * pt + r.width * s, y * pt + r.height * s)
        page.show_pdf_page(rect, src, 0)
        if letters and let:
            page.insert_text((x * pt, y * pt + 7.5), let.upper(), fontsize=8, fontname='LSB', fontfile='/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf')
    doc.save(out + '.pdf', garbage=3, deflate=True)
    pix = page.get_pixmap(dpi=dpi); pix.save(out + '.jpg', jpg_quality=92)
    page.get_pixmap(dpi=200).save(out + '_preview.png')
    return H

def svg_to_pdf(svg_text, out_pdf, w_mm, h_mm):
    from playwright.sync_api import sync_playwright
    html = f'<html><head><style>@page{{margin:0}} body{{margin:0}}</style></head><body>{svg_text}</body></html>'
    with sync_playwright() as p:
        b = p.chromium.launch(); pg = b.new_page(); pg.set_content(html)
        pg.pdf(path=out_pdf, width=f'{w_mm}mm', height=f'{h_mm}mm', print_background=True, margin=dict(top='0', bottom='0', left='0', right='0'))
        b.close()

def letter_axes(axs, xs=-0.18, y=1.06):
    """Uppercase bold panel letters A, B, C ... on a sequence of axes."""
    for i, ax in enumerate(axs): label(ax, chr(65 + i), xs, y)
