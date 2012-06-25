#! /usr/bin/env python

import os
import sys
import tempfile
import optparse
import subprocess
import re
import dendropy
from dendropy.utility import messaging

LOCAL_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
DEFAULT_SCHEMA_PATH = os.path.join(LOCAL_DIR, "schemas", "nexml", "nexml.xsd")
DEFAULT_STANDARD_VALIDATOR_PATH = os.path.join(LOCAL_DIR, "nexmlvalidator.jar")
DEFAULT_XMLLINT_PATH = "xmllint"

class XmlValidator(object):

    def __init__(self,
            schema_path=DEFAULT_SCHEMA_PATH):
        self.schema_path = schema_path
        self.error_log = []
        self.warning_log = []

    def is_valid(self):
        return len(self.error_log) == 0

    def is_invalid(self):
        return len(self.error_log) > 0

    def report_results(self, out=None):
        if out is None:
            out = sys.stderr
        for fpath, line_num, error in self.error_log:
            out.write("[%s:%s] %s\n" % (fpath, line_num, error))

class StandardNexmlValidator(XmlValidator):

    def __init__(self,
            schema_path=DEFAULT_SCHEMA_PATH):
        XmlValidator.__init__(self, schema_path=schema_path)
        pattern_template = r"\[%s\] :(\d+?):(\d+?): (.+?): (.+)"
        self.error_pattern = re.compile(pattern_template % "Error")
        self.warning_pattern = re.compile(pattern_template % "Warning")

    def validate_nexml(self, nexml_filepath):
        cmd = ["java",
                "-jar",
                DEFAULT_STANDARD_VALIDATOR_PATH,
                nexml_filepath,
                self.schema_path]
        p = subprocess.Popen(cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr.startswith("Unable to access jarfile"):
            raise OSError("Unable to access jarfile: %s" % DEFAULT_STANDARD_VALIDATOR_PATH)
        if "Failed to read schema document" in stderr:
            raise OSError("offline mode not yet supported")
        for line in stderr.split("\n"):
            m = self.error_pattern.match(line)
            if m:
                line_num, num2, cvc, msg = m.groups()
                line_num = int(line_num)
                self.error_log.append( (nexml_filepath, line_num, msg,) )
            m = self.warning_pattern.match(line)
            if m:
                line_num, num2, cvc, msg = m.groups()
                line_num = int(line_num)
                self.warning_log.append( (nexml_filepath, line_num, msg,) )

class XmllintValidator(XmlValidator):

    def __init__(self,
            schema_path=DEFAULT_SCHEMA_PATH,
            xmllint_path=DEFAULT_XMLLINT_PATH):
        XmlValidator.__init__(self, schema_path=schema_path)
        self.xmllint_path = xmllint_path
        self.error_pattern = re.compile(r"(.+?):(\d+?): (.*)")

    def validate_nexml(self, nexml_filepath):
        cmd = [self.xmllint_path,
                "--schema",
                self.schema_path,
                nexml_filepath]
        p = subprocess.Popen(cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        for line in stderr.split("\n"):
            m = self.error_pattern.match(line)
            if not m:
                continue
            fpath, line_num, error = m.groups()
            line_num = int(line_num)
            self.error_log.append( (fpath, line_num, error,) )
        try:
            d = dendropy.DataSet.get_from_path(nexml_filepath, "nexml")
        except Exception as e:
            self.error_log.append( (nexml_filepath, 0, str(e),) )


def main():
    # description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)
    # usage = "%prog [options] TREES-FILE [TREES-FILE [TREES-FILE [...]]"
    parser = optparse.OptionParser()
    parser.add_option("-s", "--schema",
            action="store",
            dest="schema",
            default=DEFAULT_SCHEMA_PATH,
            help="path to schema (default='%default')")
    parser.add_option("-x", "--xmllint",
            action="store_true",
            dest="run_xmllint",
            default=False,
            help="run xmllint as parser")
    parser.add_option("-p", "--xmllint-path",
            action="store",
            dest="xmllint_path",
            default=DEFAULT_XMLLINT_PATH,
            help="path to xmllint program (default='%default')")
    parser.add_option("-v", "--verbosity",
            action="store",
            dest="verbosity",
            type="int",
            default=3,
            help="control noise level")

    (opts, args) = parser.parse_args()
    if opts.verbosity >= 3:
        messaging_level = messaging.ConsoleMessenger.INFO_MESSAGING_LEVEL
    elif opts.verbosity >= 2:
        messaging_level = messaging.ConsoleMessenger.WARNING_MESSAGING_LEVEL
    else:
        messaging_level = messaging.ConsoleMessenger.ERROR_MESSAGING_LEVEL
    messenger = messaging.ConsoleMessenger(name="nexmlvalidator", messaging_level=messaging_level)
    if len(args) == 0:
        sys.exit("No NeXML files to validate specified")
    for arg in args:
        if opts.run_xmllint:
            validator = XmllintValidator(
                    schema_path=opts.schema,
                    xmllint_path=opts.xmllint_path)
            parser = "xmllint"
        else:
            validator = StandardNexmlValidator(schema_path=opts.schema)
            parser = "nexmlvalidator.jar"
        validator.validate_nexml(arg)
        if validator.is_invalid():
            messenger.send_error("[%s]: %s failed to validate with %d errors" % (parser, arg, len(validator.error_log)), wrap=False)
            if opts.verbosity >= 3:
                for fpath, line_num, error in validator.error_log:
                    messenger.send_info("    % 6d: %s" % (line_num, error), wrap=False)
        else:
            messenger.send_info("[%s]: %s validates with %d errors" % (parser, arg, len(validator.error_log)), wrap=False)

if __name__ == "__main__":
    main()

