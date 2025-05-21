# -*- coding: utf-8 -*-
"""
A tkinter-based UI for picking PRFs, replacing the previous PySide6 implementation.
"""

import sys
import os
import glob
import argparse
import tkinter as tk
from tkinter import filedialog, messagebox
import tkinter.font as tkFont
from os.path import exists, dirname, join

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from seispy.pickrf.pickfigure import RFFigure
from seispy.pickrf.rpickfigure import RPickFigure


class TKPickUI:
    def __init__(self, rfpath, only_r=False, xlim=[-2, 30], order='baz'):
        self.rfpath = rfpath
        self.only_r = only_r
        self.xlim = xlim
        self.order = order

        # initialize root and apply high-DPI/font scaling for high-resolution displays
        self.root = tk.Tk()
        # determine native DPI (pixels per inch)
        native_dpi = self.root.winfo_fpixels('1i')
        # compute scaling factor (relative to standard 72 DPI)
        self.scaling = native_dpi / 64.0
        # apply scaling to all Tk widgets
        self.root.tk.call('tk', 'scaling', self.scaling)

        default_size = int(10 * self.scaling) if 'scaling' in locals() else 12
        self.btn_font = tkFont.Font(family="Arial", size=default_size, weight="bold")

        self.root.title("PickRF")
        # Optionally set window icon if you have a .ico file
        icon_path = join(dirname(__file__), 'data', 'seispy.ico')
        if os.path.exists(icon_path):
            self.root.iconphoto(False, tk.PhotoImage(file=icon_path))

        self.zoom_in_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "zoom_in.png"))
        self.zoom_out_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "zoom_out.png"))
        self.back_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "back.png"))
        self.next_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "next.png"))
        self.preview_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "preview.png"))
        self.finish_icon = tk.PhotoImage(file=join(dirname(dirname(__file__)), "data", "finish.png"))

        self._build_menu()
        self._build_controls()
        self._build_canvas()
        self._bind_keys()
        self._center_window()

        # initialize figure
        self._init_figure()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_menu(self):
        menubar = tk.Menu(self.root)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Save", accelerator="Ctrl+S", command=self.plot_save)
        menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=menubar)
        # increase default font size for menus
        try:
            default_font = tkFont.nametofont("TkMenuFont")
            default_font.configure(size=int(default_font.cget("size") * self.scaling))
            tkFont.nametofont("TkDefaultFont").configure(size=int(tkFont.nametofont("TkDefaultFont").cget("size") * self.scaling))
        except Exception:
            pass

    def _build_canvas(self):
        self.fig_frame = tk.Frame(self.root)
        self.fig_frame.pack(fill=tk.BOTH, expand=1)
        # placeholder for canvas; actual fig created in _init_figure
        self.canvas = None

    def _build_controls(self):
        ctrl_frame = tk.Frame(self.root)
        ctrl_frame.pack(fill=tk.X, side=tk.TOP, pady=5)

        # amplitude controls
        amp_frame = tk.Frame(ctrl_frame)
        amp_frame.pack(side=tk.LEFT, padx=5)
        enlarge_btn = tk.Button(amp_frame, text="Amp Enlarge", command=self.enlarge,
                                image=self.zoom_in_icon, bg="#F0F0F0", font=self.btn_font,
                                relief="groove", bd=4, compound="left" )
        reduce_btn = tk.Button(amp_frame, text="Amp Reduce", command=self.reduce,
                               image=self.zoom_out_icon, bg="#F0F0F0", font=self.btn_font,
                               relief="groove", bd=4, compound="left")
        enlarge_btn.pack(side=tk.LEFT, padx=2)
        reduce_btn.pack(side=tk.LEFT, padx=2)

        # navigation controls
        nav_frame = tk.Frame(ctrl_frame)
        nav_frame.pack(side=tk.RIGHT, padx=5)
        back_btn = tk.Button(nav_frame, text="Back (z)", command=self.previous_connect,
                             image=self.back_icon, bg="#F0F0F0", font=self.btn_font,
                             relief="groove", bd=4, compound="left")
        next_btn = tk.Button(nav_frame, text="Next (c)", command=self.next_connect,
                             image=self.next_icon, bg="#F0F0F0", font=self.btn_font,
                             relief="groove", bd=4, compound="left")
        preview_btn = tk.Button(nav_frame, text="Preview (Space)", command=self.plot_ui,
                                image=self.preview_icon, bg="#F0F0F0", font=self.btn_font,
                                relief="groove", bd=4, compound="left")
        finish_btn = tk.Button(nav_frame, text="Finish", command=self.finish,
                               image=self.finish_icon, bg="#F0F0F0", font=self.btn_font,
                               relief="groove", bd=4, compound="left")
        back_btn.pack(side=tk.LEFT, padx=2)
        next_btn.pack(side=tk.LEFT, padx=2)
        preview_btn.pack(side=tk.LEFT, padx=2)
        finish_btn.pack(side=tk.LEFT, padx=2)

    def _bind_keys(self):
        self.root.bind('<z>', lambda e: self.previous_connect())
        self.root.bind('<c>', lambda e: self.next_connect())
        self.root.bind('<space>', lambda e: self.plot_ui())
        self.root.bind('<Control-s>', lambda e: self.plot_save())

    def _center_window(self):
        self.root.update_idletasks()
        w = self.root.winfo_screenwidth()
        h = self.root.winfo_screenheight()
        fw = int(w * 0.9)
        fh = int(h * 0.9)
        x = (w - fw) // 2
        y = (h - fh) // 2
        self.root.geometry(f"{fw}x{fh}+{x}+{y}")

    def _init_figure(self):
        plt.rcParams['axes.unicode_minus'] = False
        # enlarge matplotlib text for high-DPI screens
        plt.rcParams['font.size'] = 10 * self.scaling
        plt.rcParams['axes.titlesize'] = 12 * self.scaling
        plt.rcParams['axes.labelsize'] = 10 * self.scaling
        plt.rcParams['xtick.labelsize'] = 8 * self.scaling
        plt.rcParams['ytick.labelsize'] = 8 * self.scaling
        if self.only_r:
            self.rffig = RPickFigure(self.rfpath,
                                     width=21, height=11,
                                     dpi=100, xlim=self.xlim)
        else:
            self.rffig = RFFigure(self.rfpath,
                                  width=21, height=11,
                                  dpi=100, xlim=self.xlim)
        self.rffig.init_canvas(order=self.order)
        # embed in Tk
        self.canvas = FigureCanvasTkAgg(self.rffig.fig, master=self.fig_frame)
        self.canvas.draw_idle()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)
        self.canvas.mpl_connect('button_press_event', self.on_click)

    def on_click(self, event):
        self.rffig.onclick(event)
        self.canvas.draw_idle()

    def previous_connect(self):
        self.rffig.butprevious()
        self.canvas.draw_idle()

    def next_connect(self):
        self.rffig.butnext()
        self.canvas.draw_idle()

    def enlarge(self):
        self.rffig.enlarge()
        self.canvas.draw_idle()

    def reduce(self):
        self.rffig.reduce()
        self.canvas.draw_idle()

    def finish(self):
        self.rffig.finish()
        self.root.destroy()

    def plot_ui(self):
        if getattr(self.rffig, 'plotfig', None) is not None:
            plt.close(self.rffig.plotfig)
        self.rffig.plot()

    def plot_save(self):
        default_name = 'R_bazorder' if self.only_r else 'RT_bazorder'
        default_fname = os.path.join(os.getcwd(), f"{self.rffig.staname}{default_name}")
        filetypes = [("PDF Files", "*.pdf"), ("PNG Images", "*.png"), ("All Files", "*.*")]
        fname = filedialog.asksaveasfilename(
            parent=self.root,
            title="Save the figure",
            initialfile=default_fname,
            filetypes=filetypes
        )
        if not fname:
            return
        if not hasattr(self.rffig, 'plotfig') or self.rffig.plotfig is None:
            self.rffig.plot()
        try:
            self.rffig.plotfig.savefig(fname, dpi=500, bbox_inches='tight')
            self.rffig.log.RFlog.info(f"Figure saved to {fname}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))
            self.rffig.log.RFlog.error(str(e))

    def _on_close(self):
        self.root.destroy()

    def start(self):
        self.root.mainloop()


def main():
    parser = argparse.ArgumentParser(description="User interface for picking PRFs (tkinter version)")
    parser.add_argument('rf_path', type=str, help='Path to PRFs')
    parser.add_argument('-a', dest='order', default='baz', metavar='baz|dis|date',
                        help="Arrangement of RFs, defaults to 'baz'")
    parser.add_argument('-r', dest='only_r', action='store_true', 
                        help="Only plot R component")
    parser.add_argument('-x', dest='xlim', default=None, type=float,
                        help="Set x-axis max limit; defaults to 30s for RT, 85s for R.")
    args = parser.parse_args()

    rfpath = args.rf_path
    if not exists(rfpath):
        raise FileNotFoundError(f"No such directory: {rfpath}")

    # decide only_r and default xlim
    sac_t = glob.glob(join(rfpath, '*_T.sac'))
    if len(sac_t) == 0:
        only_r = True
        xlim = 85
    else:
        only_r = False
        xlim = 30
    if args.only_r:
        only_r = True
    if args.xlim is not None:
        xlim = args.xlim

    app = TKPickUI(rfpath, only_r=only_r, xlim=[-2, xlim], order=args.order)
    app.start()
    sys.exit(0)

if __name__ == '__main__':
    main()