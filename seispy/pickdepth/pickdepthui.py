# import folium
from map_type import MapChoice
from get_depth import GoodDepth
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget, QHBoxLayout,\
                              QVBoxLayout, QSizePolicy, QGroupBox,\
                              QPushButton, QTableView, QLineEdit,\
                              QLabel
from PySide6.QtGui import QStandardItemModel, QStandardItem
# from PySide6.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import io
import sys


def get_color(depth, cmap='jet_r'):
    colors = cm.get_cmap(cmap, 20)
    # print([v*255 for v in colors(depth)[0:3]])
    return '#{:02x}{:02x}{:02x}'.format(*[int(v*255) for v in colors(depth)[0:3]])


class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, stack_data_path='', idx=0,
                depmin=30, depmax=60, logger=None, width=7.2, height=11, dpi=100):
        plt.rcParams['axes.unicode_minus'] = False 
        self.gooddepth = GoodDepth(stack_data_path, width=width, height=height, dpi=dpi)
        self.gooddepth.get_dep(depmin, depmax)
        self.gooddepth.bin_idx = idx
        self.gooddepth.plot_bin()
        FigureCanvas.__init__(self, self.gooddepth.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        FigureCanvas.setCursor(self, Qt.CrossCursor)


class MapUI(QWidget):
    def __init__(self, stack_data_path, depmin=30, depmax=65,
                 width=6, height=11, dpi=100):
        super(MapUI, self).__init__()
        self.map_para = {'width': width, 'height':height, 'dpi':dpi}
        self.mpl = MyMplCanvas(self, stack_data_path=stack_data_path,
                               depmin=depmin, depmax=depmax, idx=1000,
                               **self.map_para)
        self.mpl.mpl_connect('motion_notify_event', self.on_move)
        mean_la = np.mean(self.mpl.gooddepth.ccp_data.bin_loca[:, 0])
        mean_lo = np.mean(self.mpl.gooddepth.ccp_data.bin_loca[:, 1])
        # self.map_type = MapChoice()
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.mpl, 2)
        # self.maplayout = QVBoxLayout()
        # self.layout.addLayout(self.maplayout)
        self.add_layout()
        self.setLayout(self.layout)
    
    def add_layout(self):
        self.add_page_layout()
        self.add_bin_loca_layout()
        ctrlbox = QVBoxLayout()
        ctrlbox.addStretch(1)
        ctrlbox.addWidget(self.page_box)
        ctrlbox.addWidget(self.bin_loca_box)
        self.layout.addLayout(ctrlbox)

    def add_page_layout(self):
        self.page_box = QGroupBox()
        page_layout = QHBoxLayout()
        page_layout.addStretch(1)
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

    def set_table(self):
        bin_num = self.mpl.gooddepth.ccp_data.bin_loca.shape[0]
        self.model=QStandardItemModel(bin_num, 2)
        self.model.setHorizontalHeaderLabels(['Latitude', 'Longitude'])
        for i in range(bin_num):
            itemlat = QStandardItem('{:.4f}'.format(self.mpl.gooddepth.ccp_data.bin_loca[i, 0]))
            itemlon = QStandardItem('{:.4f}'.format(self.mpl.gooddepth.ccp_data.bin_loca[i, 1]))
            self.model.setItem(i, 0, itemlat)
            self.model.setItem(i, 1, itemlon)
        self.tableWidget = QTableView()
        self.tableWidget.setModel(self.model)
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)

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
    
    def on_load(self):
        self.mpl.gooddepth.plot_bin()
        self.mpl.draw()
        self.tableWidget.selectRow(self.mpl.gooddepth.bin_idx)

    def on_idx_edit_changed(self):
        text = self.idx_edit.text()
        try:
            self.mpl.gooddepth.bin_idx = int(text)-1
        except:
            self.idx_edit.setText('')
    
    def on_move(self, event):
        self.mpl.gooddepth.on_move(event)
        self.mpl.draw()
                    
def main():
    # maptype = 'OpenTopoZhLa'
    app = QApplication(sys.argv)
    stack_data_path = '/Users/xumijian/Researches/NETibetHu/stack_data_3D_noele.npz'
    mapui = MapUI(stack_data_path)
    mapui.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()