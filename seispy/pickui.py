import sys
import random
import matplotlib
import argparse

matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QSizePolicy, QWidget, QDesktopWidget, \
                            QPushButton, QHBoxLayout, QLineEdit
from numpy import arange, sin, pi
from os.path import exists
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from seispy.pickfigure import RFFigure

class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, rfpath='', width=21, height=11, dpi=100):

        plt.rcParams['axes.unicode_minus'] = False 

        print(rfpath)
        self.rffig = RFFigure(rfpath, width=width, height=height, dpi=dpi)

        FigureCanvas.__init__(self, self.rffig.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class MatplotlibWidget(QWidget):
    def __init__(self, rfpath, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.initUi(rfpath)

    def initUi(self, rfpath):
        self.layout = QVBoxLayout(self)
        self.add_btn()
        print(rfpath)
        self.mpl = MyMplCanvas(self, rfpath=rfpath, width=21, height=11, dpi=100)
        self.layout.addWidget(self.mpl, 2)
        self.mpl.mpl_connect('button_press_event', self.on_click)

        self._set_geom_center()
    
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

    def _set_geom_center(self, height=0.8, width=0.8):
        screen_resolution = QDesktopWidget().screenGeometry()
        screen_height = screen_resolution.height()
        screen_width = screen_resolution.width()
        frame_height = int(screen_height * height)
        frame_width = int(screen_width * width)

        self.setGeometry(0, 0, frame_width, frame_height)
        self.move((screen_width / 2) - (self.frameSize().width() / 2),
                  (screen_height / 2) - (self.frameSize().height() / 2)) 

    def add_btn(self):
        pre_btn = QPushButton("Previous")
        pre_btn.clicked.connect(self.previous_connect)
        next_btn = QPushButton("Next")
        next_btn.clicked.connect(self.next_connect)
        plot_btn = QPushButton("Plot")
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
        areduce_btn.clicked.connect(self.enlarge)
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
    arg = parser.parse_args()
    rfpath = arg.rf_path
    if not exists(rfpath):
        raise('No such directory of {}'.format(rfpath))
    app = QApplication(sys.argv)
    ui = MatplotlibWidget(rfpath)
    ui.show()
    sys.exit(app.exec_()) 

if __name__ == '__main__':
    main()
