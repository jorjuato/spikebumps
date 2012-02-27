import sys as sys
import pylab as pl

# epydoc -- specify the input format.
# __docformat__ = "restructuredtext en"

def simple(funcName):
    """Create a simple plot and save the plot to a .png file
    @param funcName: The name of a function, e.g. sin, cos, tan, ...
    """
    t = pl.arange(0.0, 1.0+0.01, 0.01)
    funcStr = 'pl.%s(2*2*pl.pi*t)' % (funcName,)
    s = eval(funcStr)
    pl.plot(t, s)
    pl.xlabel('time (s)')
    pl.ylabel('voltage (mV)')
    pl.title('About as simple as it gets, folks')
    pl.grid(True)
    pl.savefig('simple_plot')
    pl.show()


def usage():
    print 'Usage: python dave_simple_plot.py <func_name>'
    print 'Examples:'
    print '    python dave_simple_plot.py sin'
    print '    python dave_simple_plot.py cos'
    print '    python dave_simple_plot.py tan'
    sys.exit(-1)


def main():
    args = sys.argv[1:]
    if len(args) != 1:
        usage()
    simple(args[0])


if __name__ == '__main__':
    main()
