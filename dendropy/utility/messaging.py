#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Messaging, logging and support.
"""

import sys
import os
import logging
import textwrap
from dendropy.utility import metavar

###############################################################################
## metavar.LOGGING

def get_logging_level():
    if metavar.LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG
        elif os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[metavar.LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET
    return level

def get_logger(name="dendropy"):
    """
    Returns a logger with name set as given, and configured
    to the level given by the environment variable metavar.LOGGING_LEVEL_ENVAR.
    """

#     package_dir = os.path.dirname(module_path)
#     config_filepath = os.path.join(package_dir, metavar.LOGGING_CONFIG_FILE)
#     if os.path.exists(config_filepath):
#         try:
#             logging.config.fileConfig(config_filepath)
#             logger_set = True
#         except:
#             logger_set = False
    logger = logging.getLogger(name)
    if not hasattr(logger, 'is_configured'):
        logger.is_configured = False
    if not logger.is_configured:
        level = get_logging_level()
        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        simple_formatter = logging.Formatter("%(levelname) 8s: %(message)s")
        raw_formatter = logging.Formatter("%(message)s")
        default_formatter = None
        logging_formatter = default_formatter
        if metavar.LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[metavar.LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[metavar.LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[metavar.LOGGING_FORMAT_ENVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
        logger.is_configured = True
    return logger

class ConsoleMessenger(object):

    ERROR_MESSAGING_LEVEL = 2000
    WARNING_MESSAGING_LEVEL = 1000
    INFO_MESSAGING_LEVEL = 100

    def __init__(self,
            name="DendroPy",
            messaging_level=None,
            dest=sys.stderr):
        self.name = name
        if messaging_level is None:
            self.messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
        else:
            self.messaging_level = messaging_level
        self.primary_out = dest
        self.text_wrapper = textwrap.TextWrapper(width=78, initial_indent= "", subsequent_indent=" " * (len(self.name) + 2))
        self.message_leader = {
                ConsoleMessenger.ERROR_MESSAGING_LEVEL : self.error_leader,
                ConsoleMessenger.WARNING_MESSAGING_LEVEL : self.warning_leader,
                ConsoleMessenger.INFO_MESSAGING_LEVEL : self.info_leader
                }
        self.silent = False

    def error_leader(self):
        return self.name + ": [ERROR] "

    def warning_leader(self):
        return self.name + ": [WARNING] "

    def info_leader(self):
        return self.name + ": "

    def format_message(self, msg, level, wrap=True, prefix=""):
        if not wrap:
            msg = self.message_leader[level]() + prefix + msg
        else:
            if prefix:
                full_leader = self.message_leader[level]() + prefix
                msg = textwrap.fill(
                    msg,
                    width=self.text_wrapper.width,
                    initial_indent=full_leader,
                    subsequent_indent=" " * len(full_leader))
            else:
                msg = self.text_wrapper.fill(self.message_leader[level]() + msg)
        return msg

    def log(self, msg, level=0, wrap=True, prefix="", newline=True):
        if self.silent:
            return
        if level >= self.messaging_level:
            msg = self.format_message(msg, level, wrap=wrap, prefix=prefix)
            self.primary_out.write(msg)
            if newline:
                self.primary_out.write("\n")

    def log_lines(self, msg, level=None, wrap=True, prefix=""):
        for line in msg:
            self.log(msg=line, level=level, wrap=wrap, prefix=prefix)

    def error(self, msg, wrap=True, prefix=""):
        self.log(msg, level=ConsoleMessenger.ERROR_MESSAGING_LEVEL, wrap=wrap, prefix=prefix)

    def warning(self, msg, wrap=True, prefix=""):
        self.log(msg, level=ConsoleMessenger.WARNING_MESSAGING_LEVEL, wrap=wrap, prefix=prefix)

    def info(self, msg, wrap=True, prefix=""):
        self.log(msg, level=ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap, prefix=prefix)

    def info_lines(self, msg, wrap=True, prefix=""):
        for line in msg:
            self.log(msg=line, level=ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap, prefix=prefix)

    def info_raw(self, msg):
        if self.messaging_level <= ConsoleMessenger.INFO_MESSAGING_LEVEL:
            self.primary_out.write(self.info_leader() + msg)




