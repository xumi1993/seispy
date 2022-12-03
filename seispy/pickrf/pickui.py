import sys
import os
import argparse
# matplotlib.use("Qt5Agg")
from PySide6.QtGui import QIcon, QKeySequence, QAction, QGuiApplication, QShortcut
from PySide6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, \
                            QSizePolicy, QWidget, \
                            QPushButton, QHBoxLayout, QFileDialog            
from os.path import exists, dirname, join
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from seispy.pickrf.pickfigure import RFFigure
from seispy.pickrf.rpickfigure import RPickFigure
import glob


class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, rfpath='', only_r=False, width=21, height=11,
                 dpi=100, xlim=[-2, 30], order='baz'):

        plt.rcParams['axes.unicode_minus'] = False 

        if only_r:
            self.rffig = RPickFigure(rfpath, width=width, height=height, dpi=dpi, xlim=xlim)
            self.rffig.init_canvas(order=order)
        else:
            self.rffig = RFFigure(rfpath, width=width, height=height, dpi=dpi, xlim=xlim)
            self.rffig.init_canvas(order=order)

        FigureCanvas.__init__(self, self.rffig.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Policy.Expanding,
                                   QSizePolicy.Policy.Expanding)
        FigureCanvas.updateGeometry(self)


class MatplotlibWidget(QMainWindow):
    def __init__(self, rfpath, only_r=False, xlim=[-2, 30], 
                 order='baz', parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.only_r = only_r
        self.initUi(rfpath, only_r, xlim, order=order)

    def initUi(self, rfpath, only_r, xlim, order='baz'):
        self.layout = QVBoxLayout()
        self.add_btn()
        self.mpl = MyMplCanvas(self, rfpath=rfpath, only_r=only_r, width=21, height=11,
                               dpi=100, xlim=xlim, order=order)
        self.layout.addWidget(self.mpl, 2)
        self.mpl.mpl_connect('button_press_event', self.on_click)

        main_frame = QWidget()
        self.setCentralWidget(main_frame)
        main_frame.setLayout(self.layout)

        saveAction = QAction('&Save', self)        
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save this figure')
        saveAction.triggered.connect(self.plot_save)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(saveAction)

        self._set_geom_center()
        self._define_global_shortcuts()
        self.setWindowTitle('PickRF')
        self.setWindowIcon(QIcon(join(dirname(__file__), 'data', 'seispy.png')))

    def on_click(self, event):
        self.mpl.rffig.onclick(event)
        self.mpl.draw()

    def previous_connect(self):
        self.mpl.rffig.butprevious()
        self.mpl.draw()

    def next_connect(self):
        self.mpl.rffig.butnext()
        self.mpl.draw()

    def enlarge(self):
        self.mpl.rffig.enlarge()
        self.mpl.draw()

    def reduce(self):
        self.mpl.rffig.reduce()
        self.mpl.draw()

    def finish(self):
        self.mpl.rffig.finish()
        QApplication.quit()

    def plot_ui(self):
        if self.mpl.rffig.plotfig is not None:
            plt.close(self.mpl.rffig.plotfig)
        self.mpl.rffig.plot()

    def plot_save(self):
        if self.only_r:
            default_name = 'R_bazorder'
        else:
            default_name = 'RT_bazorder'
        fileName_choose, filetype = QFileDialog.getSaveFileName(self,
                                    "Save the figure",
                                    os.path.join(os.getcwd(), self.mpl.rffig.staname + default_name), 
                                    "PDF Files (*.pdf);;Images (*.png);;All Files (*)")

        if fileName_choose == "":
            return
        if not hasattr(self.mpl.rffig, 'plotfig'):
            self.mpl.rffig.plot()
        try:
            self.mpl.rffig.plotfig.savefig(fileName_choose, dpi=500, bbox_inches='tight')
            self.mpl.rffig.log.RFlog.info('Figure saved to {}'.format(fileName_choose))
        except Exception as e:
            self.mpl.rffig.log.RFlog.error('{}'.format(e))

    def _set_geom_center(self, height=1, width=1):
        screen_resolution = QGuiApplication.primaryScreen().geometry()
        screen_height = screen_resolution.height()
        screen_width = screen_resolution.width()
        frame_height = int(screen_height * height)
        frame_width = int(screen_width * width)

        self.setGeometry(0, 0, frame_width, frame_height)
        self.move(int((screen_width / 2) - (self.frameSize().width() / 2)),
                  int((screen_height / 2) - (self.frameSize().height() / 2)))

    def _define_global_shortcuts(self):
        self.key_c = QShortcut(QKeySequence('c'), self)
        self.key_c.activated.connect(self.next_connect)
        self.key_z = QShortcut(QKeySequence('z'), self)
        self.key_z.activated.connect(self.previous_connect)
        self.key_space = QShortcut(QKeySequence('Space'), self)
        self.key_space.activated.connect(self.plot_ui)

    def add_btn(self):
        pre_btn = QPushButton("Back (z)")
        pre_btn.clicked.connect(self.previous_connect)
        next_btn = QPushButton("Next (c)")
        next_btn.clicked.connect(self.next_connect)
        plot_btn = QPushButton("Preview (Space)")
        plot_btn.clicked.connect(self.plot_ui)
        finish_btn = QPushButton("Finish")
        finish_btn.clicked.connect(self.finish)
        btnbox = QHBoxLayout()
        btnbox.addStretch(1)
        btnbox.addWidget(pre_btn)
        btnbox.addWidget(next_btn)
        btnbox.addWidget(plot_btn)
        btnbox.addWidget(finish_btn)

        enlarge_btn = QPushButton("Amp enlarge")
        enlarge_btn.clicked.connect(self.enlarge)
        areduce_btn = QPushButton("Amp reduce")
        areduce_btn.clicked.connect(self.reduce)
        pathbox = QHBoxLayout()
        pathbox.addWidget(enlarge_btn)
        pathbox.addWidget(areduce_btn)

        ctrl_layout = QHBoxLayout()
        ctrl_layout.addLayout(pathbox)
        ctrl_layout.addLayout(btnbox)

        self.layout.addLayout(ctrl_layout)


def main():
    parser = argparse.ArgumentParser(description="User interface for picking PRFs")
    parser.add_argument('rf_path', type=str, help='Path to PRFs')
    parser.add_argument('-a', help='Arrangement of RFs, defaults to \'baz\'', dest='order',
                        default='baz', type=str, metavar='baz|dis')
    parser.add_argument('-x', help="Set x limits of the current axes, defaults to 30s for RT, 85s for R.",
                        dest='xlim', default=None, type=float, metavar='xmax')
    arg = parser.parse_args()
    rfpath = arg.rf_path
    if not exists(rfpath):
        raise FileNotFoundError('No such directory of {}'.format(rfpath))
    if len(glob.glob(join(rfpath, '*_T.sac'))) == 0:
        only_r = True
        xlim = 85
    else:
        only_r = False
        xlim = 30
    if arg.xlim is not None:
        xlim = arg.xlim
    app = QApplication(sys.argv)
    ui = MatplotlibWidget(rfpath, only_r=only_r, xlim=[-2, xlim], order=arg.order)
    ui.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
