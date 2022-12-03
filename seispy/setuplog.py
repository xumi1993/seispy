import logging
from os.path import expanduser, join


class setuplog(object):
    def __init__(self, filename=join(expanduser('~'), '.RF.log')):
        self.filename = filename
        fh = logging.FileHandler(filename)
        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        self.RF2depthlog = logging.getLogger('RF2depth')
        if not self.RF2depthlog.handlers:
            self.RF2depthlog.setLevel(logging.INFO)
            self.RF2depthlog.addHandler(ch)
        self.RFlog = logging.getLogger('RF')
        if not self.RFlog.handlers:
            self.RFlog.setLevel(logging.INFO)
            # self.RFlog.removeHandler(ch)
            self.RFlog.addHandler(ch)
        self.Batlog = logging.getLogger('Bat')
        if not self.Batlog.handlers:
            self.Batlog.setLevel(logging.INFO)
            # self.Batlog.removeHandler(ch)
            # self.Batlog.removeHandler(fh)
            self.Batlog.addHandler(ch)
            self.Batlog.addHandler(fh)
        self.CCPlog = logging.getLogger('CCP')
        if not self.CCPlog.handlers:
            self.CCPlog.setLevel(logging.INFO)
            self.CCPlog.addHandler(ch)
        self.ModCreatorlog = logging.getLogger('ModCreator')
        if not self.ModCreatorlog.handlers:
            self.ModCreatorlog.setLevel(logging.INFO)
            self.ModCreatorlog.addHandler(ch)
        self.PickDepthlog = logging.getLogger('PickDepth')
        if not self.PickDepthlog.handlers:
            self.PickDepthlog.setLevel(logging.INFO)
            # self.PickDepthlog.addHandler(ch)


if __name__ == '__main__':
    logger = setuplog()
    logger.RFlog.info('print')