#!/usr/bin/python
import numpy
from numpy import *
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def plot_points(trn_points, dev_points, file):
    plt.figure()
    x = trn_points[:, 0]
    y = trn_points[:, 1]
    plt.plot(x, y, 'r+', label='train')
    x1 = dev_points[:, 0]
    y1 = dev_points[:, 1]
    plt.plot(x1, y1, 'bo', label='dev')
    plt.legend(loc='best')
    plt.title('The scattering of points')
    plt.xlabel('First dimension')
    plt.ylabel('Second dimension')
    plt.savefig(file, format='eps')

def plot_graph(points):
    (min_x, min_y) = numpy.min(points, axis=0)
    (max_x, max_y) = numpy.max(points, axis=0)
    print 'range of first dimension is ( %.3f, %.3f )' % (min_x, max_x)
    print 'range of second dimension is ( %.3f, %.3f )' % (min_y, max_y)
    (floor_x, floor_y, ceil_x, ceil_y) = (floor(min_x), floor(min_y), ceil(max_x), ceil(max_y))
    x_interval = (ceil_x - floor_x) / 5
    y_interval = (ceil_y - floor_y) / 5
    x_grid = arange(floor_x, ceil_x, x_interval)
    y_grid = arange(floor_y, ceil_y, y_interval)
    X, Y = meshgrid(x_grid, y_grid)
    plt.figure()
    plt.plot(points[:,0], points[:,1], 'r+')
    plt.show()

def plot_likelihood(likelihoods, file, title):
    plt.figure()
    x = array(range(0, len(likelihoods)))
    x *= 10
    #likelihoods /= 1000000
    y = likelihoods
    plt.plot(x, y, 'r')
    plt.title(title)
    plt.xlim((0, len(likelihoods) * 10))
    margin_size = 0.1 * (numpy.max(likelihoods) - numpy.min(likelihoods))
    plt.ylim((-margin_size + (numpy.min(likelihoods)), numpy.max(likelihoods) + margin_size))
    plt.xlabel('Iteration #(10)')
    plt.ylabel('BLEU-score')
    plt.savefig(file, format='eps')

def plot_comparison_likelihood(diff_likelihoods, labels, file, title):
    plt.figure()
    n_curves = len(diff_likelihoods)
    n_iteration = len(diff_likelihoods[0])

    for i in xrange(n_curves):
        likelihoods = diff_likelihoods[i]
        assert len(likelihoods) == n_iteration, 'Run for different iterations'
        curr_label = labels[i]
        x = array(range(0, len(likelihoods)))
        likelihoods /= 10000
        y = likelihoods
        #plt.plot(x, y, label=curr_label)
        plt.plot(x, y, '--', label=curr_label)
    plt.legend(loc='best')
    plt.title(title)
    plt.xlim((0, n_iteration+1))
    margin_size = 0.1 * (numpy.max(diff_likelihoods) - numpy.min(diff_likelihoods))
    plt.ylim((-margin_size + (numpy.min(diff_likelihoods)), numpy.max(diff_likelihoods) + margin_size))
    plt.xlabel('Iteration #')
    plt.ylabel('Log likelihood')
    plt.savefig(file, format='eps')

def draw_likehood_from_file(filename, output_file):
    f = open(filename)
    likelihoods = []
    line = f.readline()
    while line is not None:
        line = line.strip()
        if line == '':
            break
        score = float(line)
        likelihoods.append(score)
        line = f.readline()
    plot_likelihood(array(likelihoods), output_file, 'iter-bleu-score')

if __name__ == '__main__':
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    #output_file = sys.argv[3]
    f1 = open(file1)
    f2 = open(file2)
    likelihoods1 = []
    likelihoods2 = []
    #line1 = f1.readline()
    #line2 = f2.readline()
    #line1 = f1.readline()
    #line2 = f2.readline()
    #line3 = f3.readline()
    #line1 = f1.readline()
    #line2 = f2.readline()
    #line3 = f3.readline()
    for i in xrange(0, 161):
        line1 = f1.readline().strip()
        line2 = f2.readline().strip()
        if line1 == '':
            break
        score1 = float(line1)
        score2 = float(line2)
        likelihoods1.append(score1)
        likelihoods2.append(score2)
    labels=[]
    likelihoods=[]
    likelihoods.append(array(likelihoods1[3:]))
    likelihoods.append(array(likelihoods2[3:]))
    labels.append('anneal')
    labels.append('non-anneal')
    plot_comparison_likelihood(likelihoods, labels, 'simulated-annealing.eps', 'SHRG sampling with or without simulated annealing')
    #draw_likehood_from_file(sys.argv[1], sys.argv[2])
    #plot_likelihood(array(likelihoods), output_file, 'iter-likelihood')
