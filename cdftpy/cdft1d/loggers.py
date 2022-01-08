import io
import logging

simple_format = logging.Formatter('%(message)s')


def get_stream_handler():
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(simple_format)
    return stream_handler

def get_stream_logger(name, capture=True, quiet = False, debug=False):

    logger = logging.getLogger(name)
    logger.addHandler(logging.NullHandler())
    if not quiet:
        console_handler = get_stream_handler()
        logger.addHandler(console_handler)

    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if capture:
        stream_handler = get_stream_handler()
        streamer = io.StringIO()
        stream_handler.setStream(streamer)
        logger.addHandler(stream_handler)
    else:
        streamer = None

    logger.propagate = False

    return logger, streamer

