def main():

    import matplotlib
    matplotlib.use('Qt4Agg')
    import numpy
    import pylab

    t_list, x_cm_list, y_cm_list = numpy.loadtxt('cm_history.txt',unpack=True)
    pylab.plot(t_list, y_cm_list)
    pylab.show()

if __name__ == '__main__':

    main()
