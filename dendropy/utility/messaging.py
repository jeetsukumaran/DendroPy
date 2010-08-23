#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Messaging, logging and support.
"""

import sys
import os
import logging
import textwrap

###############################################################################
## LOGGING

_LOGGING_LEVEL_ENVAR="DENDROPY_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR="DENDROPY_LOGGING_FORMAT"

def get_logging_level():
    if _LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET
    return level

def get_logger(name="dendropy"):
    """
    Returns a logger with name set as given, and configured
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """

#     package_dir = os.path.dirname(module_path)
#     config_filepath = os.path.join(package_dir, _LOGGING_CONFIG_FILE)
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
        if _LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "NONE":
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

def deprecation(message, logger_obj=None, stacklevel=3):
    try:
        import warnings
        warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)
    except:
        if logger_obj:
            logger_obj.warning(message)

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
        self.text_wrapper = textwrap.TextWrapper(width=78, subsequent_indent=" " * (len(self.name) + 2))
        self.message_leader = {
                ConsoleMessenger.ERROR_MESSAGING_LEVEL : self.error_leader,
                ConsoleMessenger.WARNING_MESSAGING_LEVEL : self.warning_leader,
                ConsoleMessenger.INFO_MESSAGING_LEVEL : self.info_leader
                }

    def error_leader(self):
        return self.name + ": <<< ERROR >>> "

    def warning_leader(self):
        return self.name + ": [[[ WARNING ]]] "

    def info_leader(self):
        return self.name + ": "

    def format_message(self, msg, level, wrap=True):
        msg = self.message_leader[level]() + msg
        if wrap:
            msg = self.text_wrapper.fill(msg)
        return msg

    def send(self, msg, level, wrap=True, newline=True):
        if level is None:
            level = ConsoleMessenger.INFO_MESSAGING_LEVEL
        if level >= self.messaging_level:
            msg = self.format_message(msg, level, wrap=wrap)
            self.primary_out.write(msg)
            if newline:
                self.primary_out.write("\n")

    def send_lines(self, msg, level=None, wrap=True, prefix=""):
        for line in msg:
            self.send(msg=prefix+line, level=level, wrap=wrap)

    def send_error(self, msg, wrap=True):
        self.send(msg, level=ConsoleMessenger.ERROR_MESSAGING_LEVEL, wrap=wrap)

    def send_warning(self, msg, wrap=True):
        self.send(msg, level=ConsoleMessenger.WARNING_MESSAGING_LEVEL, wrap=wrap)

    def send_info(self, msg, wrap=True):
        self.send(msg, level=ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap)

    def send_info_lines(self, msg, wrap=True, prefix=""):
        for line in msg:
            self.send(msg=prefix+line, level=ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap)

    def write(self, msg):
        self.primary_out.write(self.info_leader() + msg)
