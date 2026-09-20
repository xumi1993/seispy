import logging
from os.path import expanduser, join

_DEFAULT_LOG_FILE = join(expanduser('~'), '.RF.log')


class SetupLog:
    default_level = logging.INFO
    default_logs={
        "Invlog":("Inv",default_level,"stream_handler"),
        "RF2depthlog":("RF2depth",default_level,"stream_handler"),
        "RFlog":("RF",default_level,"stream_handler"),
        "Batlog":("Bat",default_level,"stream_handler","file_handler"),
        "CCPlog":("CCP",default_level,"stream_handler"),
        "ModCreatorlog":("ModCreator",default_level,"stream_handler"),
        "PickDepthlog": ("PickDepth", default_level,"stream_handler")
    }
    def __init__(self, filename=_DEFAULT_LOG_FILE):
        """Create configured SeisPy loggers.

        :param filename: Log file for branches using a file handler. None
            initializes stream handlers without opening or truncating a file
        :type filename: str or None, optional
        """
        self.filename = filename
        fh = None if filename is None else logging.FileHandler(filename, mode='w')
        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s')
        if fh is not None:
            fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        for loger_branch, config in self.default_logs.items():
            # init, setlevel
            log = logging.getLogger(config[0])
            log.setLevel(config[1])

            # add handler
            if not log.hasHandlers():
                if "file_handler" in config and fh is not None:
                    log.addHandler(fh)
                if "stream_handler" in config:
                    log.addHandler(ch)

            # attach to class
            setattr(self,loger_branch,log)

def inversion_logger(log=None):
    """Return the SeisPy inversion logger without reopening the RF log file.

    :param log: Existing SetupLog, logging.Logger or LoggerAdapter; defaults to None
    :type log: object, optional
    :return: Logger used for inversion progress and diagnostics
    :rtype: logging.Logger or logging.LoggerAdapter
    """
    if log is None:
        logger = logging.getLogger(SetupLog.default_logs['Invlog'][0])
        if logger.hasHandlers():
            return logger
        return SetupLog(filename=None).Invlog
    if isinstance(log, SetupLog):
        return log.Invlog
    if isinstance(log, (logging.Logger, logging.LoggerAdapter)):
        return log
    raise TypeError('log must be a SetupLog, logging.Logger or LoggerAdapter')


if __name__ == '__main__':
    logger = SetupLog()
    logger.RFlog.info('print')
