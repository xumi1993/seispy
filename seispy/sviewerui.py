import sys
import os
import argparse
# matplotlib.use("Qt5Agg")
from PyQt5.QtGui import QIcon, QKeySequence
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, \
                            QSizePolicy, QWidget, QDesktopWidget, \
                            QPushButton, QHBoxLayout, QFileDialog, \
                            QAction, QShortcut
from os.path import exists, dirname, join
from seispy.setuplog import setuplog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt



class SFigure(Figure):
    def __init__(self, eqs, para, width=21, height=11, dpi=100):
        self.eqs = eqs
        self.row_num = eqs.shape[0]
        self.para = para
        self.log = setuplog()
        self.idx = 0
        self.init_figure(width=width, height=height, dpi=dpi)
        self.plot()
        self.drop_lst = []
        self.drop_color = 'lightgray'

    def init_figure(self, width=21, height=11, dpi=100):
        self.fig, self.axs = plt.subplots(3, 1, sharex=True, figsize=(width, height), dpi=dpi, tight_layout=True)
    
    def set_properties(self):
        date = self.eqs.iloc[self.idx]['date'].strftime("%Y.%m.%dT%H:%M:%S")
        dis = self.eqs.iloc[self.idx]['dis']
        bazi = self.eqs.iloc[self.idx]['bazi']
        mag = self.eqs.iloc[self.idx]['mag']
        title = '({}/{}) {}, Baz: {:.2f}$^{{\circ}}$, Distance:{:.2f} km, Mw: {:.1f}'.format(self.idx+1, self.row_num, date, bazi, dis, mag)
        self.fig.suptitle(title)
        for ax in self.axs:
            ax.set(xlim=[-self.para.time_before, self.para.time_after])
            ax.minorticks_on()
        self.axs[2].set_xlabel('Time (s)')

    def clear(self):
        for ax in self.axs:
            ax.cla()

    def plot(self, lc='tab:blue'):
        st = self.eqs.iloc[self.idx]['data'].rf
        self.time_axis = st[0].times()-self.para.time_before
        for i in range(3):
            self.axs[i].plot(self.time_axis, st[i].data, color=lc)
            self.axs[i].axvline(x=0, lw=1, ls='--', color='r')
            self.axs[i].set_ylabel(st[i].stats.channel)
        self.set_properties()

    def next_action(self):
        self.clear()
        self.idx += 1
        if self.idx >= self.row_num:
            self.idx = 0
        if self.eqs.index[self.idx] in self.drop_lst:
            self.plot(lc=self.drop_color)
        else:
            self.plot()

    def back_action(self):
        self.clear()
        self.idx -= 1
        if self.idx < 0 :
            self.idx = 0
        if self.eqs.index[self.idx] in self.drop_lst:
            self.plot(lc=self.drop_color)
        else:
            self.plot()

    def drop(self):
        if self.eqs.index[self.idx] in self.drop_lst:
            return
        self.drop_lst.append(self.eqs.index[self.idx])
        self.plot(lc=self.drop_color)

    def cancel(self):
        if self.eqs.index[self.idx] not in self.drop_lst:
            return
        self.drop_lst.remove(self.eqs.index[self.idx])
        self.plot()

    def finish(self):
        self.eqs.drop(self.drop_lst, inplace=True)

class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, eqs=None, para=None, width=21, height=11, dpi=100):

        plt.rcParams['axes.unicode_minus'] = False 
        self.sfig = SFigure(eqs, para, width=width, height=height, dpi=dpi)

        FigureCanvas.__init__(self, self.sfig.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class MatplotlibWidget(QMainWindow):
    def __init__(self, eqs, para, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.initUi(eqs, para)

    def initUi(self, eqs, para):
        self.mpl = MyMplCanvas(self, eqs=eqs, para=para, width=21, height=11, dpi=100)
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.mpl, 2)

        main_frame = QWidget()
        self.setCentralWidget(main_frame)
        main_frame.setLayout(self.layout)

        self._set_geom_center()
        self._define_global_shortcuts()
        self.setWindowTitle('PickSPhease')
        self.setWindowIcon(QIcon(join(dirname(__file__), 'data', 'seispy.png')))

    def exit_app(self):
        self.close()

    def next_connect(self):
        self.mpl.sfig.next_action()
        self.mpl.draw()
    
    def back_connect(self):
        self.mpl.sfig.back_action()
        self.mpl.draw()
    
    def drop_connect(self):
        self.mpl.sfig.drop()
        self.mpl.draw()
    
    def cancel_connect(self):
        self.mpl.sfig.cancel()
        self.mpl.draw()
    
    def finish_connect(self):
        self.mpl.sfig.finish()
        QApplication.quit()
    
    def _define_global_shortcuts(self):
        self.key_c = QShortcut(QKeySequence('c'), self)
        self.key_c.activated.connect(self.next_connect)
        self.key_z = QShortcut(QKeySequence('z'), self)
        self.key_z.activated.connect(self.back_connect)
        self.key_d = QShortcut(QKeySequence('d'), self)
        self.key_d.activated.connect(self.drop_connect)
        self.key_a = QShortcut(QKeySequence('a'), self)
        self.key_a.activated.connect(self.cancel_connect)
        self.key_enter = QShortcut(QKeySequence('Return'), self)
        self.key_enter.activated.connect(self.finish_connect)

    def _set_geom_center(self, height=0.7, width=1):
        screen_resolution = QDesktopWidget().screenGeometry()
        screen_height = screen_resolution.height()
        screen_width = screen_resolution.width()
        frame_height = int(screen_height * height)
        frame_width = int(screen_width * width)

        self.setGeometry(0, 0, frame_width, frame_height)
        self.move((screen_width / 2) - (self.frameSize().width() / 2),
                  (screen_height / 2) - (self.frameSize().height() / 2))