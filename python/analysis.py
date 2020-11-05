import numpy as np
import matplotlib.pyplot as plt

from util import load_results, load_dist, load_Tlimits

eps = 1e-10

def compare_test_results(dataset, model_name):
    print('Comparing test results for {} on {} graphs:'.format(model_name, dataset))
    s1, t1 = load_results('../data/results/{}_{}.csv'.format(dataset, 'default'))
    s2, t2 = load_results('../data/results/{}_{}.csv'.format(dataset, model_name))

    s1, t1, s2, t2 = remove_uncompleted(s1, t1, s2, t2)
    s1 = np.array(s1)
    t1 = np.array(t1)
    s2 = np.array(s2)
    t2 = np.array(t2)
    print('Time difference: {:.4f}'.format(np.sum(t1-t2)))
    print('Steps difference: {}'.format(np.sum(s1-s2)))

    print('Time ratio: {:.4f}'.format(np.sum(t1)/np.sum(t2)))
    print('Steps ratio: {:.4f}'.format(np.sum(s1)/np.sum(s2)))

    print('Avg time ratio: {:.4f}'.format(np.mean(t1/t2)))
    print('Avg steps ratio: {:.4f}'.format(np.mean(s1/s2)))
    print()

def remove_uncompleted(s1, t1, s2, t2):
    n_removed = 0
    i = 0
    while i != len(s1):
        if s1[i] < 0 or t1[i] < 0 or s2[i] < 0 or t2[i] < 0:
            s1.pop(i)
            s2.pop(i)
            t1.pop(i)
            t2.pop(i)
            n_removed += 1
            continue
        i += 1
    print("Number of uncompleted:", n_removed)
    return s1, t1, s2, t2

def plot_distribution(x, plot_name='Tlimit_dist'):
    #x += eps
    # histogram on linear scale
    plt.subplot(211)
    hist, bins, _ = plt.hist(x, bins=8)
    plt.ylabel('N')
    plt.xlabel('Tlimit')

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.subplot(212)
    plt.hist(x, bins=logbins)
    plt.xscale('log')
    plt.ylabel('N')
    plt.xlabel('Tlimit')
    #plt.show()
    plt.savefig('figures/' + plot_name + '.png')

def plot_time_dist(path):
    Tlimits, times = load_dist(path)
    # filter out failed Tlimits
    Tlimits = [Tlimits[i] for i in range(len(times)) if times[i] > 0]
    times = [times[i] for i in range(len(times)) if times[i] > 0]
    
    plt.figure(figsize=(10,5))
    plt.subplot(211)
    plt.ylabel('time (s)')
    plt.xlabel('Tlimit')
    plt.plot(Tlimits, times)
    plt.axvline(x=0.025, color='r')

    plt.subplot(212)
    plt.xscale('log')
    plt.ylabel('time (s)')
    plt.xlabel('Tlimit')
    plt.plot(Tlimits, times)
    plt.axvline(x=0.025, color='r')

    plt.savefig('figures/' + plot_name + '.png')
    #plt.show()

if __name__ == "__main__":
    
    
    plot_distribution(load_Tlimits('../datasets/docking_train_tlimits.csv'), plot_name='Tlimit_dist_docking')
    plot_distribution(load_Tlimits('../datasets/product_train_tlimits.csv'), plot_name='Tlimit_dist_product')
    plot_distribution(load_Tlimits('../datasets/rand_train_tlimits.csv'), plot_name='Tlimit_dist_rand')
    
    



    #plot_time_dist()

    #compare_test_results("protein", "gnn")
    #compare_test_results("dense", "gnn")
    #compare_test_results("rand", "gnn")
    #compare_test_results("dimacs", "gnn")

    #compare_test_results("rand", "gbr")
    #compare_test_results("rand", "gnn")
    #compare_test_results("rand", "gin")
    #compare_test_results("rand", "adaboost")

    #compare_test_results("dimacs", "gbr")
    #compare_test_results("dimacs", "gnn")
    #compare_test_results("dimacs", "gin")
    #compare_test_results("dimacs", "adaboost")

    #_, Tlimits = load_data('../data/train_data_3.csv')
    #print('Number of examples with invalid Tlimits:', sum([t <= 0 for t in Tlimits])
    #plot_distribution(Tlimits)

