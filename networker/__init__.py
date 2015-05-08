__version__ = "0.1.0"

# initialize logging
import logging 

logger = logging.getLogger('networker')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# need to set this, otherwise 'root' logger also logs
logging.getLogger('networker').propagate = False
