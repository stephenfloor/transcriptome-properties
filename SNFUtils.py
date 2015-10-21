#!/usr/bin/env python 

import os,subprocess,sys
from itertools import izip
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr, spearmanr 

def plot_txs (tx_i, tx_j, clust_i, clust_j, 
              x_i, y_i, error_i, xfit_i, yfit_i,
              x_j, y_j, error_j, xfit_j, yfit_j):
    
    #plt.figure()
    fig, ax = plt.subplots()
    plt.errorbar(x_i, y_i, marker='x', yerr=error_i, ls="None")
    plt.plot(xfit_i, yfit_i)
    plt.errorbar(x_j, y_j, marker='s', yerr=error_j, ls="None")
    plt.plot(xfit_j, yfit_j)
    
    diffs = [yfit_i[i] - yfit_j[i] for i in range(len(yfit_i))]
    plt.plot(xfit_i, diffs)

    ymax = max( max(y_i), max(yfit_i), max(y_j), max(yfit_j), max(diffs))*1.2+1
    ymin = min( min(y_i), min(yfit_i), min(y_j), min(yfit_j), min(diffs))*1.2-1
    #plt.plot(x, y, 'x', xnew, ynew)
    plt.axis([.5, 9.5, ymin, ymax])
    ax.set_xticks(range(1,10))
    ax.set_xticklabels(['80S', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto']) 
    plt.legend([tx_i + "_clust %s" % clust_i, tx_i + "_fit", tx_j + "_clust %s" % clust_j, tx_j + "_fit", "difference"], fontsize=8)

    #residuals = sum(infodict['fvec']**2)
    #plt.title("Sq. Resid.: %5.4f; Res/Median: %5.4f" % (residuals, residuals/ymedian))
    #perr = np.sqrt(np.diag(pcov))
    #perr_percent = [ np.fabs(perr[i]/popt[i]) for i in range(len(popt))]
    #avg_percent_error = np.mean(perr_percent)
    #total_percent_error = sum(perr_percent)
    #weighted_perr = sum([ perr_percent[i] * np.fabs(popt[i]) for i in range(len(popt))])
    
    prsn = pearsonr( yfit_i, yfit_j)[0]
    sprmn = spearmanr( yfit_i, yfit_j)[0]

    plt.text(.75, 1, "Parms pearson: %3.2f spearman: %3.2f" % (prsn, sprmn), fontsize=8)
    plt.savefig("%s vs %s" % (tx_i, tx_j))
    plt.close(fig) 
              
def third_order_poly_fit_plot (x, y, outname, yerror):
    def func(x, p1, p2, p3, p4):
        return p1 + p2 * x + p3 * x**2 + p4 * x**3
    
    xdata = np.array(x) 
    ydata = np.array(y)
    ymedian = np.median(y) 
    xnew = np.arange(1, max(x), 0.001) 
    popt, pcov, infodict, mesg, ier = curve_fit(func, xdata, ydata,p0=(1, 1, 1, 1),full_output=1) 
    ynew = [func(i, popt[0], popt[1], popt[2], popt[3]) for i in xnew]
    #plt.figure()
    fig, ax = plt.subplots()
    plt.errorbar(x, y, marker='x', yerr=yerror, ls="None")
    plt.plot(xnew, ynew)
    #plt.plot(x, y, 'x', xnew, ynew)
    plt.axis([.5, 9.5, 0, max( max(y), max(ynew) ) + 1])
    ax.set_xticklabels(['', '80S', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto']) 
    plt.legend(['Input', 'Third order polynomial'])
    residuals = sum(infodict['fvec']**2)
    plt.title("Sq. Resid.: %5.4f; Res/Median: %5.4f" % (residuals, residuals/ymedian))
    perr = np.sqrt(np.diag(pcov))
    perr_percent = [ np.fabs(perr[i]/popt[i]) for i in range(len(popt))]
    avg_percent_error = np.mean(perr_percent)
    total_percent_error = sum(perr_percent)
    weighted_perr = sum([ perr_percent[i] * np.fabs(popt[i]) for i in range(len(popt))])
    prsn = pearsonr( [func(i, popt[0], popt[1], popt[2], popt[3]) for i in x], y)[0]
    sprmn = spearmanr( [func(i, popt[0], popt[1], popt[2], popt[3]) for i in x], y)[0]
    def prt(inp): #"pretty" 
        return ["%3.2f" % inp[i] for i in range(len(inp))]

    plt.text(.75, 1, "Parms %s\nerrors %s\n%% error %s\nmean %%: %3.2f sum %%: %3.2f weighted %%: %3.2f pearson: %3.2f spearman: %3.2f" % 
              (prt(popt), prt(perr), prt(perr_percent), avg_percent_error, total_percent_error, weighted_perr, prsn, sprmn),
              fontsize=8)
    plt.savefig(outname)
    plt.close(fig) 

def plot_dist (x, y, outname, yerror):
    
    xdata = np.array(x) 
    ydata = np.array(y)
    ymedian = np.median(y) 

    fig, ax = plt.subplots()
    plt.errorbar(x, y, marker='x', yerr=yerror, ls="None")

    plt.axis([.5, 9.5, 0, max(y) + 1])
    ax.set_xticks(range(1,10))
    ax.set_xticklabels(['80S', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto']) 
    plt.legend(['Input'])

    plt.savefig(outname)
    plt.close(fig) 

def plot_dist_fancy (x, y, outname, yerror, title):
    
    xdata = np.array(x) 
    ydata = np.array(y)
    ymedian = np.median(y) 


    fig, ax = plt.subplots()
    plt.errorbar(x, y, marker='o', markersize=16, color='k', yerr=yerror, ls="None")



    # plot them all
    #plt.axis([.5, len(x) + 0.5, 0, max(max(y), max([y[i] + yerror[i] for i in range(len(y))]))])
    #ax.set_xticks(range(1,len(x)+1))
    #ax.set_xticklabels(['40S', '60S', '80S', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto'], size=20, rotation=45) 

    # leave out 40/60 for comparison to clustering 
    plt.axis([2.5, len(x) + 0.5, 0, max(max(y[2:]), max([y[i] + yerror[i] for i in range(2, len(y))]))])
    ax.set_xticks(range(3,len(x)+1))
    ax.set_xticklabels(['80S', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto'], size=20, rotation=45) 

    # leave out 40/60/80 for comparison to frac-seq rtpcr 
    #plt.axis([3.5, len(x) + 0.5, 0, max(max(y), max([y[i] + yerror[i] for i in range(len(y))]))])
    #ax.set_xticks(range(4,len(x)+1))
    #ax.set_xticklabels(['poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'cyto'], size=20, rotation=45) 

    ylabels = ax.get_yticks().tolist() 
    ax.set_yticklabels(ylabels, size=20)

    plt.title(title, size=24)
    plt.ylabel("TPM", size=24)
    
    plt.tick_params(which='both', length=8, width=2, pad=10)

    plt.tick_params(
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',
    right='off')         # ticks along the top edge are off
    #labelbottom='off') # labels along the bottom edge are off


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['bottom'].set_linewidth(4)
    ax.spines['left'].set_linewidth(4)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)

    plt.axvline(x=10.5, linewidth=5, color='#cccccc', dashes=(4,16), dash_capstyle="round") 

    plt.tight_layout()

    plt.savefig(outname)
    plt.close(fig) 


def pairwise(t):
    it = iter(t)
    return izip(it, it)

def chunkwise(t, size=2):
    it = iter(t)
    return izip(*[it]*size)

def stdout_from_command(command):
    p = subprocess.Popen(command,
                         stdout = subprocess.PIPE,
                         shell = True)
    return iter(p.stdout.readline, b'')

def safe_open_file(filename):
    if (os.path.exists(filename)):
        sys.exit("FATAL: file %s exists; cowardly refusing to overwrite." % filename)
    try:
        outfile = open(filename, "w")
    except:
        sys.exit("FATAL: cannot open file %s for writing." % filename)

    return outfile

def prompt(promptstr): 
    print promptstr

    inp = raw_input("\nContinue?  [y/n]").lower()

    if (inp == "n" or inp == "no"):
        sys.exit(0)

def is_number(s):
    if s is None: 
        return False
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False

def merge(d1, d2, merge_fn=lambda x,y:y):
    """
    http://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression
    Merges two dictionaries, non-destructively, combining 
    values on duplicate keys as defined by the optional merge
    function.  The default behavior replaces the values in d1
    with corresponding values in d2.  (There is no other generally
    applicable merge strategy, but often you'll have homogeneous 
    types in your dicts, so specifying a merge technique can be 
    valuable.)

    Examples:

    >>> d1
    {'a': 1, 'c': 3, 'b': 2}
    >>> merge(d1, d1)
    {'a': 1, 'c': 3, 'b': 2}
    >>> merge(d1, d1, lambda x,y: x+y)
    {'a': 2, 'c': 6, 'b': 4}

    """
    print "---- d1 ----" 
    print d1
    print "---- d2 ----" 
    print d2 
    result = dict(d1)
    for k,v in d2.iteritems():
        if k in result:
            result[k] = merge_fn(result[k], v)
        else:
            result[k] = v
    return result
