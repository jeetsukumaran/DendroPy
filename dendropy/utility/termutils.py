import sys

def ttysize():
    try:
        fp = os.popen('stty -a', 'r')
        ln1 = fp.readline()
        fp.close()
        if not ln1:
            raise ValueError, 'tty size not supported for input'
        vals = {'rows':None, 'columns':None}
        for ph in string.split(ln1, ';'):
            x = string.split(ph)
            if len(x) == 2:
                vals[x[0]] = x[1]
                vals[x[1]] = x[0]
        return vals['rows'], vals['columns']
    except:
        return 40, 80

def posix_terminal_width():
    """Return estimated terminal width."""
    width = 0
    try:
        import struct, fcntl, termios
        s = struct.pack('HHHH', 0, 0, 0, 0)
        x = fcntl.ioctl(1, termios.TIOCGWINSZ, s)
        width = struct.unpack('HHHH', x)[1]
    except IOError:
        pass
    if width <= 0:
        try:
            width = int(os.environ['COLUMNS'])
        except:
            pass
    if width <= 0:
        width = 80
    return width

def win_terminal_width():
    from ctypes import windll, create_string_buffer
    # stdin handle is -10, stdout handle is -11, stderr handle is -12
    h = windll.kernel32.GetStdHandle(-12)
    csbi = create_string_buffer(22)
    res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
    if res:
        import struct
        (bufx, bufy, curx, cury, wattr,
         left, top, right, bottom, maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        sizex = right - left + 1
        sizey = bottom - top + 1
    else:
        sizex, sizey = 80, 25 # can't determine actual size - return default values
    return sizex-1

def terminal_width():
    if sys.platform.startswith('win'):
        return win_terminal_width()
    else:
        return posix_terminal_width()
