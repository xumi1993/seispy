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
        self.RFlog = logging.getLogger('RF')
        self.RFlog.setLevel(logging.INFO)
        self.RFlog.removeHandler(ch)
        self.RFlog.addHandler(ch)
        self.Batlog = logging.getLogger('Bat')
        self.Batlog.setLevel(logging.INFO)
        self.Batlog.removeHandler(ch)
        self.Batlog.removeHandler(fh)
        self.Batlog.addHandler(ch)
        self.Batlog.addHandler(fh)


if __name__ == '__main__':
    logger = setuplog()
    logger.RFlog.info('print')