#!/usr/bin/env
import tempfile
import subprocess
import os

DENDROPY_USER_DIR = os.path.expanduser("~/.dendropy/")
TMP_FILE_LIST_FILENAME = os.path.join(DENDROPY_USER_DIR, "tmpfiles")
def db_show_tree(t, prefix="tmp"):
    newick = str(t)
    pd = os.path.join(DENDROPY_USER_DIR, "tmp")
    if not os.path.exists(pd):
        os.makedirs(pd)
    t = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', prefix=prefix, delete=False, dir=pd)
    n = t.name
    t.write(newick)
    t.close()
    subprocess.Popen(["open", "-a", "/Applications/FigTree v1.2.1.app/", n])
    tflf = open(TMP_FILE_LIST_FILENAME, "a")
    tflf.write("%s\n" % n)
    tflf.close()
    
def db_clean_tmp_files():
    try:
        tflf = open(TMP_FILE_LIST_FILENAME, "r")
    except:
        return
    undeleted = []
    for line in tflf:
        fn = line.strip()
        if os.path.exists(fn):
            d = os.path.dirname(fn)
            try:
                os.remove(fn)
            except:
                sys.stderr.write("Could not remove %s\n" % fn)
                undeleted.append(fn)
        if not os.path.exists(fn):
            try:
                os.removedirs(d)
            except:
                pass
    os.remove(TMP_FILE_LIST_FILENAME)
    if undeleted:
        tflf = open(TMP_FILE_LIST_FILENAME, "a")
        for n in undeleted:
            tflf.write("%s\n" % n)
        tflf.close()
    
    
    
