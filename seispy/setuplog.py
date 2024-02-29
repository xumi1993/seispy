import logging
from os.path import expanduser, join


class SetupLog(object):
    default_level = logging.INFO
    default_logs={
        "RF2depthlog":("RF2depth",default_level,"stream_handler"),
        "RFlog":("RF",default_level,"stream_handler"),
        "Batlog":("Bat",default_level,"stream_handler","file_handler"),
        "CCPlog":("CCP",default_level,"stream_handler"),
        "ModCreatorlog":("ModCreator",default_level,"stream_handler"),
        "PickDepthlog": ("PickDepth", default_level,"stream_handler")
    }
    def __init__(self, filename=join(expanduser('~'), '.RF.log')):
        """
        use default_logs to gen loggers
        change default_logs for future changes if needed,
        check logger level with logger.getEffectiveLevel

        """
        self.filename = filename
        fh = logging.FileHandler(filename, mode='w')
        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        for loger_branch, config in self.default_logs.items():
            # init, setlevel
            log = logging.getLogger(config[0])
            log.setLevel(config[1])

            # add handler
            if not log.hasHandlers():
                if "file_handler" in config:
                    log.addHandler(fh)
                if "stream_handler" in config:
                    log.addHandler(ch)

            # attach to class
            setattr(self,loger_branch,log)


if __name__ == '__main__':
    logger = SetupLog()
    logger.RFlog.info('print')