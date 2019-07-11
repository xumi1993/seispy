import sys
import os
import random
import matplotlib
import argparse

matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QSizePolicy, QWidget, QDesktopWidget, \
                            QPushButton, QHBoxLayout, QLineEdit, QFileDialog, QAction
from numpy import arange, sin, pi
from os.path import exists
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from seispy.pickfigure import RFFigure


class PlotCanvas(FigureCanvas):
    def __init__(self, parent, fig):

        plt.rcParams['axes.unicode_minus'] = False 

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class  PreviewWidget(QWidget):
    def __init__(self, fig, parent=None):
        super(PreviewWidget, self).__init__(parent)
        self.fig = fig
        self.initUi()

    def initUi(self):
        self.layout = QVBoxLayout(self)
        self.mpl = PlotCanvas(self, self.fig)
        self.layout.addWidget(self.mpl)

        self.setWindowTitle('PRFs order by back-azimuth')    
        


def show_plotui(fig):
    app = QApplication(sys.argv)
    ui = PreviewWidget(fig)
    ui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    h = plt.figure(figsize=(11.7, 8.3))
    show_plotui(h)
