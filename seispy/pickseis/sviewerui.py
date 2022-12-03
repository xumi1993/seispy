# matplotlib.use("Qt5Agg")
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QKeySequence, QGuiApplication, QShortcut, QAction
from PySide6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, \
                            QSizePolicy, QWidget, \
                            QPushButton, QHBoxLayout, QFileDialog
from os.path import exists, dirname, join
from seispy.setuplog import setuplog
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class SFigure(Figure):
    def __init__(self, eqs, para, logger, width=21, height=11, dpi=100):
        self.eqs = eqs
        self.logger = logger
        self.row_num = eqs.shape[0]
        self.picker_time = pd.DataFrame({'trigger_shift': np.array([eqs.iloc[i]['data'].trigger_shift for i in range(self.row_num)])})
        self.picker_time.set_index(eqs.index, inplace=True)
        self.para = para
        self.log = setuplog()
        self.idx = 0
        self.init_figure(width=width, height=height, dpi=dpi)
        self.plot()
        self.drop_lst = []
        self.drop_color = 'lightgray'

    def init_figure(self, width=21, height=11, dpi=100):
        self.fig, self.axs = plt.subplots(3, 1, sharex=True, figsize=(width, height), dpi=dpi, tight_layout=True)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.cursor[0].get_visible() != visible
        # self.horizontal_line.set_visible(visible)
        for cursor in self.cursor:
            cursor.set_visible(visible)
        # self.text.set_visible(visible)
        return need_redraw
    
    def set_properties(self):
        date = self.eqs.iloc[self.idx]['date'].strftime("%Y.%m.%dT%H:%M:%S")
        lat = self.eqs.iloc[self.idx]['evla']
        lon = self.eqs.iloc[self.idx]['evlo']
        dep = self.eqs.iloc[self.idx]['evdp']
        dis = self.eqs.iloc[self.idx]['dis']
        bazi = self.eqs.iloc[self.idx]['bazi']
        mag = self.eqs.iloc[self.idx]['mag']
        title = '({}/{}) {}, Lat: {:.2f}, Lon: {:.2f}, Depth: {:.1f}, \nBaz: {:.2f}$^{{\circ}}$, Distance:{:.2f} km, Mw: {:.1f}'.format(
                self.idx+1, self.row_num, date, lat, lon, dep, bazi, dis, mag)
        self.fig.suptitle(title)
        for ax in self.axs:
            ax.set(xlim=[-self.para.time_before, self.para.time_after])
            ax.autoscale(enable=True,axis='y', tight=True)
            ax.minorticks_on()
        self.axs[2].set_xlabel('Time (s)')
        self.axs[0].legend(loc=2)


    def clear(self):
        for ax in self.axs:
            ax.cla()

    def _plot_line(self, st, lc='tab:blue'):
        self.picker = []
        self.lines = []
        for i in range(3):
            self.lines.append(self.axs[i].plot(self.time_axis, st[i].data, color=lc))
            self.axs[i].axvline(x=0, lw=1.5, color='r', label='Theoretical')
            self.picker.append(self.axs[i].axvline(self.picker_time.iloc[self.idx]['trigger_shift'],
                                                   label='Triggered', color='k', lw=1.5, visible=True))
            self.axs[i].set_ylabel(st[i].stats.channel)

    def plot(self, lc='tab:blue'):
        st = self.eqs.iloc[self.idx]['data'].st_pick
        if len(st) == 0:
            return
        # t1, t2 = self.eqs.iloc[self.idx]['data'].t1_pick, st = self.eqs.iloc[self.idx]['data'].t2_pick
        self.time_axis = st[0].times()-self.para.time_before
        self._plot_line(st, lc)
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
    
    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                for ax in self.axs:
                    ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x = event.xdata
            # update the line positions
            # self.horizontal_line.set_ydata(y)
            for cursor in self.cursor:
                cursor.set_xdata(x)
                # ax.figure.canvas.draw()
            # self.cursor[1].set_xdata(x)
            # self.axs[1].figure.canvas.draw()
        
    def on_click(self, event):
        if not event.inaxes:
            return
        if self.eqs.index[self.idx] in self.drop_lst:
            return
        self.picker_time.iloc[self.idx]['trigger_shift'] = event.xdata
        for picker in self.picker:
            picker.set_xdata(self.picker_time.iloc[self.idx]['trigger_shift'])
            picker.set_visible(True)

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
        self.logger.RFlog.info('{} is marked'.format(self.eqs.loc[self.eqs.index[self.idx], 'datestr']))
        # self.clear()
        for line in self.lines:
            line[0].set_color(self.drop_color)
        # self._plot_line(self.eqs.iloc[self.idx]['data'].st_pick, )

    def cancel(self):
        if self.eqs.index[self.idx] not in self.drop_lst:
            return
        self.drop_lst.remove(self.eqs.index[self.idx])
        self.logger.RFlog.info('{} is unmarked'.format(self.eqs.loc[self.eqs.index[self.idx], 'datestr']))
        # self.clear()
        for line in self.lines:
            line[0].set_color('tab:blue')
        # self._plot_line(self.eqs.iloc[self.idx]['data'].st_pick)

    def finish(self):
        self.eqs.drop(self.drop_lst, inplace=True)
        for i, row in self.eqs.iterrows():
            row['data'].trigger_shift = self.picker_time['trigger_shift'][i]

class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, eqs=None, para=None, logger=None, width=21, height=11, dpi=100):
        plt.rcParams['axes.unicode_minus'] = False 
        self.sfig = SFigure(eqs, para, logger, width=width, height=height, dpi=dpi)
        FigureCanvas.__init__(self, self.sfig.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Policy.Expanding,
                                   QSizePolicy.Policy.Expanding)
        FigureCanvas.updateGeometry(self)
        FigureCanvas.setCursor(self, Qt.CursorShape.CrossCursor)

class MatplotlibWidget(QMainWindow):
    def __init__(self, eqs, para, logger, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.initUi(eqs, para, logger)

    def initUi(self, eqs, para, logger):
        self.mpl = MyMplCanvas(self, eqs=eqs, para=para, logger=logger, width=8, height=4, dpi=50)
        # self.mpl.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.mpl.mpl_connect('button_press_event', self.on_click)
        self.layout = QVBoxLayout()
        self.add_btn()
        self.layout.addWidget(self.mpl, 2)

        main_frame = QWidget()
        self.setCentralWidget(main_frame)
        main_frame.setLayout(self.layout)

        self._set_geom_center()
        self._define_global_shortcuts()
        self.setWindowTitle('Pick S Phase')
        self.setWindowIcon(QIcon(join(dirname(__file__), 'data', 'seispy.png')))

    def add_btn(self):
        pre_btn = QPushButton("Back (z)")
        # pre_btn.setCursor(Qt.ArrowCursor)
        pre_btn.clicked.connect(self.back_connect)
        next_btn = QPushButton("Next (c)")
        # next_btn.setCursor(Qt.ArrowCursor)
        next_btn.clicked.connect(self.next_connect)
        finish_btn = QPushButton("Finish (Enter)")
        # finish_btn.setCursor(Qt.ArrowCursor)
        finish_btn.clicked.connect(self.finish_connect)
        btnbox = QHBoxLayout()
        btnbox.addStretch(1)
        btnbox.addWidget(pre_btn)
        btnbox.addWidget(next_btn)
        btnbox.addWidget(finish_btn)

        self.mark_btn = QPushButton("Mark to poor (d)")
        # mark_btn.setCursor(Qt.ArrowCursor)
        self.mark_btn.setCheckable(True)
        self.mark_btn.clicked.connect(self.drop_connect)
        pathbox = QHBoxLayout()
        pathbox.addWidget(self.mark_btn)

        ctrl_layout = QHBoxLayout()
        ctrl_layout.addLayout(pathbox)
        ctrl_layout.addLayout(btnbox)

        self.layout.addLayout(ctrl_layout)

    # def mark_btn_connect(self, mark_btn):
    #     if mark_btn.isChecked():
    #         self.drop_connect()
    #     else:
    #         self.cancel_connect()

    def exit_app(self):
        self.close()

    def on_click(self, event):
        self.mpl.sfig.on_click(event)
        self.mpl.draw()
    
    def on_mouse_move(self, event):
        self.mpl.sfig.on_mouse_move(event)
        self.mpl.draw()

    def check_mark_btn(self):
        sfig = self.mpl.sfig
        if sfig.eqs.index[sfig.idx] in sfig.drop_lst:
            self.mark_btn.setChecked(True)
        else:
            self.mark_btn.setChecked(False)

    def check_drop(self):
        sfig = self.mpl.sfig
        if sfig.eqs.index[sfig.idx] in sfig.drop_lst:
            self.mpl.sfig.cancel()
            self.mark_btn.setChecked(False)
        else:
            self.mpl.sfig.drop()
            self.mark_btn.setChecked(True)

    def next_connect(self):
        self.mpl.sfig.next_action()
        self.mpl.draw()
        self.check_mark_btn()
    
    def back_connect(self):
        self.mpl.sfig.back_action()
        self.mpl.draw()
        self.check_mark_btn()
    
    def drop_connect(self):
        # self.mark_btn.setChecked(True)
        # self.mpl.sfig.drop()
        self.check_drop()
        self.mpl.draw()
    
    # def cancel_connect(self):
    #     self.mark_btn.setChecked(False)
    #     self.mpl.sfig.cancel()
    #     self.mpl.draw()
    
    def finish_connect(self):
        self.mpl.sfig.finish()
        self.close()
        # QApplication.quit()
    
    def _define_global_shortcuts(self):
        self.key_c = QShortcut(QKeySequence('c'), self)
        self.key_c.activated.connect(self.next_connect)
        self.key_z = QShortcut(QKeySequence('z'), self)
        self.key_z.activated.connect(self.back_connect)
        self.key_d = QShortcut(QKeySequence('d'), self)
        self.key_d.activated.connect(self.drop_connect)
        # self.key_a = QShortcut(QKeySequence('a'), self)
        # self.key_a.activated.connect(self.cancel_connect)
        self.key_enter = QShortcut(QKeySequence('Return'), self)
        self.key_enter.activated.connect(self.finish_connect)

    def _set_geom_center(self, height=0.7, width=1):
        screen_resolution = QGuiApplication.primaryScreen().geometry()
        screen_height = screen_resolution.height()
        screen_width = screen_resolution.width()
        frame_height = int(screen_height * height)
        frame_width = int(screen_width * width)

        self.setGeometry(0, 0, frame_width, frame_height)
        self.move((screen_width / 2) - (self.frameSize().width() / 2),
                  (screen_height / 2) - (self.frameSize().height() / 2))