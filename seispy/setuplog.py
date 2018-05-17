import logging


class setuplog(object):
    def __init__(self):
        # self.filename = filename
        # fh = logging.FileHandler(filename)
        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s')
        # fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        self.RFlog = logging.getLogger('RF')
        self.RFlog.setLevel(logging.INFO)
        self.RFlog.addHandler(ch)


if __name__ == '__main__':
    logger = setuplog()
    logger.RFlog.info('print')