# import folium
from .get_depth import GoodDepth
from seispy.setuplog import setuplog
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget, QHBoxLayout,\
                              QVBoxLayout, QSizePolicy, QGroupBox,\
                              QPushButton, QTableWidget, QLineEdit,\
                              QLabel, QAbstractItemView, QTableWidgetItem,\
                              QToolButton, QMessageBox, QFileDialog,\
                              QPlainTextEdit, QHeaderView
from PySide6.QtGui import QIcon, QShortcut, QKeySequence, QGuiApplication
# from PySide6.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys
from os.path import dirname, join, abspath, isdir, isfile
import os
import logging
import argparse


class LoggerText(logging.Handler):
    def __init__(self, parent):
        super(LoggerText, self).__init__()
        self.widget = parent
        self.widget.setReadOnly(True)

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)


def get_color(depth, cmap='jet_r'):
    colors = cm.get_cmap(cmap, 20)
    # print([v*255 for v in colors(depth)[0:3]])
    return '#{:02x}{:02x}{:02x}'.format(*[int(v*255) for v in colors(depth)[0:3]])


class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, stack_data_path='', idx=0,
                depmin=30, depmax=60, logger=None, width=7.2,
                height=11, dpi=100, smooth=10):
        plt.rcParams['axes.unicode_minus'] = False 
        self.gooddepth = GoodDepth(stack_data_path, logger, width=width,
                                   height=height, dpi=dpi, smooth=smooth)
        self.gooddepth.get_dep(depmin, depmax)
        self.gooddepth.bin_idx = idx
        self.gooddepth._get_next_bin()
        self.gooddepth.plot_bin()
        FigureCanvas.__init__(self, self.gooddepth.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # FigureCanvas.setCursor(self, Qt.CrossCursor)


class MapUI(QWidget):
    def __init__(self, stack_data_path, depmin=30, depmax=65,
                 width=6, height=11, dpi=100, idx=0, smooth=10):
        super(MapUI, self).__init__()
        self.log = setuplog()
        # self._set_geom_center()
        self.add_log_layout()
        self.map_para = {'width': width, 'height':height, 'dpi':dpi}
        self.mpl = MyMplCanvas(self, stack_data_path=stack_data_path,
                               depmin=depmin, depmax=depmax, idx=idx,
                               logger=self.log, smooth=smooth, **self.map_para)
        self.mpl.mpl_connect('button_press_event', self.on_click)
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.mpl, 7)
        self.add_layout()
        self.setLayout(self.layout)
        self._define_global_shortcuts()
        self.setWindowTitle('Pick Depth')
        self.setWindowIcon(QIcon(join(dirname(dirname(abspath(__file__))), 'data', 'seispy.png')))
    
    def _set_geom_center(self, height=1, width=1):
        screen_resolution = QGuiApplication.primaryScreen().geometry()
        screen_height = screen_resolution.height()
        screen_width = screen_resolution.width()
        frame_height = int(screen_height * height)
        frame_width = int(screen_width * width)

        self.setGeometry(0, 0, frame_width, frame_height)
        self.move(int((screen_width / 2) - (self.frameSize().width() / 2)),
                  int((screen_height / 2) - (self.frameSize().height() / 2)))

    def add_layout(self):
        self.add_save_layout()
        self.add_page_layout()
        self.add_bin_loca_layout()
        # self.add_log_layout()
        ctrlbox = QVBoxLayout()
        # ctrlbox.addStretch(1)
        ctrlbox.addWidget(self.save_box)
        ctrlbox.addWidget(self.page_box)
        ctrlbox.addWidget(self.bin_loca_box)
        ctrlbox.addWidget(self.log_box)
        self.layout.addLayout(ctrlbox, 3)

    def add_save_layout(self):
        self.save_box = QGroupBox('Database saving')
        save_layout = QHBoxLayout()
        self.openbutton = QToolButton()
        self.openbutton.setIcon(QIcon(join(dirname(dirname(abspath(__file__))), 'data', 'saveicon.png')))
        self.openbutton.clicked.connect(self.on_save_file)
        self.save_edit = QLineEdit()
        self.save_edit.setText(join(os.getcwd(), 'good_depths.dat'))
        self.save_edit.setFocusPolicy(Qt.ClickFocus)
        self.savebutton = QPushButton('Save')
        self.savebutton.clicked.connect(self.on_save)
        save_layout.addWidget(self.openbutton, 0)
        save_layout.addWidget(self.save_edit, 9)
        save_layout.addWidget(self.savebutton, 3)
        self.save_box.setLayout(save_layout)

    def add_page_layout(self):
        self.page_box = QGroupBox('Page up/down')
        page_layout = QHBoxLayout()
        # page_layout.addStretch(1)
        self.puButton = QPushButton()
        self.puButton.setText('Back (z)')
        self.puButton.clicked.connect(self.page_up)
        self.pdButton = QPushButton()
        self.pdButton.setText('Next (c)')
        self.pdButton.clicked.connect(self.page_down)
        page_layout.addWidget(self.puButton)
        page_layout.addWidget(self.pdButton)
        self.page_box.setLayout(page_layout)
    
    def add_bin_loca_layout(self):
        self.bin_loca_box = QGroupBox('Bin locations')
        layout = QVBoxLayout()
        show_layout = QHBoxLayout()
        label_bar = QLabel()
        label_bar.setText('Index: ')
        self.idx_edit = QLineEdit()
        self.idx_edit.setText('{}'.format(self.mpl.gooddepth.bin_idx+1))
        self.idx_edit.textChanged.connect(self.on_idx_edit_changed)
        self.idx_edit.setFocusPolicy(Qt.ClickFocus)
        self.loadButton = QPushButton()
        self.loadButton.setText('Load')
        self.loadButton.clicked.connect(self.on_load)
        show_layout.addWidget(label_bar)
        show_layout.addWidget(self.idx_edit)
        show_layout.addWidget(self.loadButton)
        self.set_table()
        layout.addLayout(show_layout)
        layout.addWidget(self.tableWidget)
        self.bin_loca_box.setLayout(layout)

    def add_log_layout(self):
        self.log_box = QGroupBox('Logs')
        layout = QHBoxLayout()
        self.log_text = QPlainTextEdit()
        loghandler = LoggerText(self.log_text)
        self.log.PickDepthlog.addHandler(loghandler)
        # self.log_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.log_text)
        self.log_box.setLayout(layout)

    def set_table(self):
        bin_num = self.mpl.gooddepth.ccp_data.bin_loca.shape[0]
        self.tableWidget = QTableWidget(bin_num, 3)
        self.tableWidget.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Depth'])
        header = self.tableWidget.horizontalHeader()  
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)
        self.tableWidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        for i in range(bin_num):
            itemlat = QTableWidgetItem('{:.3f}'.format(self.mpl.gooddepth.ccp_data.bin_loca[i, 0]))
            itemlon = QTableWidgetItem('{:.3f}'.format(self.mpl.gooddepth.ccp_data.bin_loca[i, 1]))
            itemdep = QTableWidgetItem('{:.2f}'.format(self.mpl.gooddepth.good_depth.iloc[i]['depth']))
            self.tableWidget.setItem(i, 0, itemlat)
            self.tableWidget.setItem(i, 1, itemlon)
            self.tableWidget.setItem(i, 2, itemdep)
        self.tableWidget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tableWidget.cellPressed.connect(self.on_select_row)

    def page_up(self):
        self.mpl.gooddepth.page_up()
        self.mpl.draw_idle()
        self.idx_edit.setText('{}'.format(self.mpl.gooddepth.bin_idx+1))
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)

    
    def page_down(self):
        self.mpl.gooddepth.page_down()
        self.mpl.draw_idle()
        self.idx_edit.setText('{}'.format(self.mpl.gooddepth.bin_idx+1))
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)
    
    def on_select_row(self, event):
        idx = self.tableWidget.currentIndex().row()
        self.mpl.gooddepth.bin_idx = idx
        self.idx_edit.setText('{}'.format(idx+1))

    def on_load(self):
        self.mpl.gooddepth.plot_bin()
        self.mpl.draw()
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)

    def on_save_file(self):
        fedit = self.save_edit.text()
        fname = QFileDialog.getSaveFileName(self,
                'Save as',
                dirname(fedit),
                'All Files (*);; ASCII File (*.dat)')
        self.save_edit.setText(fname[0])

    def on_save(self):
        fname = self.save_edit.text()
        if isdir(fname):
            msg = QMessageBox.critical(self, 'Error',
                 'Cannot create file. {} is a directory.'.format(fname),
                 QMessageBox.Ok)
            return
        if isfile(fname):
            msg = QMessageBox(self)
            msg.setWindowTitle('Warning')
            msg.setText('{} already exists. Do you want to overwrite it.'.format(fname))
            msg.setStandardButtons(QMessageBox.StandardButton.No|QMessageBox.StandardButton.Yes)
            msg.setIcon(QMessageBox.Icon.Warning)
            button = msg.exec()
            if button == QMessageBox.StandardButton.No:
                return
        try:
            self.mpl.gooddepth.write(self.save_edit.text())
            msg = QMessageBox.information(self, 'Information',
                'Successfully save database to {}'.format(fname),
                QMessageBox.Ok)
        except Exception as e:
            msg = QMessageBox.critical(self, 'Error',
                'Error: {}.'.format(e),
                QMessageBox.Ok)
            return

    def on_idx_edit_changed(self):
        text = self.idx_edit.text()
        try:
            self.mpl.gooddepth.bin_idx = int(text)-1
        except:
            self.idx_edit.setText('')
    
    def on_click(self, event):
        bin_idx = self.mpl.gooddepth.bin_idx
        self.mpl.gooddepth.on_click(event)
        itemdep = QTableWidgetItem('{:.2f}'.format(self.mpl.gooddepth.good_depth.iloc[bin_idx]['depth']))
        self.tableWidget.setItem(bin_idx, 2, itemdep)
        self.mpl.draw()
    
    def _define_global_shortcuts(self):
        self.key_c = QShortcut(QKeySequence('c'), self)
        self.key_c.activated.connect(self.page_down)
        self.key_z = QShortcut(QKeySequence('z'), self)
        self.key_z.activated.connect(self.page_up)
        self.key_enter = QShortcut(QKeySequence('Enter'), self)
        self.key_enter.activated.connect(self.on_load)
        self.key_save = QShortcut(QKeySequence.StandardKey.Save, self)
        self.key_save.activated.connect(self.on_save)

           
def main():
    parser = argparse.ArgumentParser(description="User interface for picking PRFs")
    parser.add_argument('stack_data_path', type=str, help='Path to CCP stacked data')
    parser.add_argument('-d', help='Depth range contain target interface', metavar='dep_min/dep_max')
    parser.add_argument('-i', help='Specify starting index of bins', type=int, default=1, metavar='index')
    parser.add_argument('-s', help='Smoothing scale in km', type=float, default=2, metavar='smooth')
    args = parser.parse_args()
    try:
        dep_range = [float(v) for v in args.d.split('/')]
    except:
        raise ValueError('Error format of depth range')
    if args.i < 1:
        raise ValueError('The index should be greater than 0')
    app = QApplication(sys.argv)    
    mapui = MapUI(args.stack_data_path, depmin=dep_range[0], depmax=dep_range[1], idx=args.i-1)
    mapui.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()