import numpy as np
import matplotlib.pyplot as plt

from util import load_results, load_dist, load_Tlimits, load_predictions

eps = 1e-10

def compare_test_results(dataset, model_name):
    print('Comparing test results for {} on {} graphs:'.format(model_name, dataset))
    s1, t1 = load_results('../results/{}_{}.csv'.format(dataset, 'default'))
    s2, t2 = load_results('../results/{}_{}.csv'.format(dataset, model_name))

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
    # histogram on linear scale
    plt.subplot(211)
    hist, bins, _ = plt.hist(x, bins=8)
    plt.ylabel('N')
    #plt.xlabel('Tlimit')

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.subplot(212)
    plt.hist(x, bins=logbins)
    plt.xscale('log')
    plt.ylabel('N')
    plt.xlabel('Tlimit')
    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()

def plot_time_dist(path, plot_name):
    Tlimits, times = load_dist(path)
    # filter out failed Tlimits
    Tlimits = [Tlimits[i] for i in range(len(times)) if times[i] > 0]
    times = [times[i] for i in range(len(times)) if times[i] > 0]
    
    plt.figure(figsize=(10,5))
    plt.subplot(211)
    plt.ylabel('time (s)')
    plt.plot(Tlimits, times)
    plt.axvline(x=0.025, color='r')

    plt.subplot(212)
    plt.xscale('log')
    plt.ylabel('time (s)')
    plt.xlabel('Tlimit')
    plt.plot(Tlimits, times)
    plt.axvline(x=0.025, color='r')

    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()
    #plt.show()

def plot_time_dist_models(dataset, path, plot_name, ix):
    Tlimits, times = load_dist(path)
    # filter out failed Tlimits
    Tlimits = [Tlimits[i] for i in range(len(times)) if times[i] > 0]
    times = [times[i] for i in range(len(times)) if times[i] > 0]
    
    plt.figure(figsize=(10,5))
    #plt.subplot(211)
    #plt.ylabel('time (s)')
    #plt.plot(Tlimits, times)
    #plt.axvline(x=0.025, color='r')
    
    #plt.subplot(212)
    plt.xscale('log')
    plt.ylabel('time (s)')
    plt.xlabel('Tlimit')
    plt.plot(Tlimits, times, label='_nolegend_')


    # get predictions for each model
    p_xgb = load_predictions('../pred/{}_xgb.csv'.format(dataset))[ix]
    p_gcn = load_predictions('../pred/{}_gcn.csv'.format(dataset))[ix]
    p_gan = load_predictions('../pred/{}_gat.csv'.format(dataset))[ix]
    p_gin = load_predictions('../pred/{}_gin.csv'.format(dataset))[ix]
    p_svr = load_predictions('../pred/{}_svr_wl.csv'.format(dataset))[ix]

    plt.axvline(x=0.025, color='r') # default
    plt.axvline(x=p_xgb, color='y') # XGB
    plt.axvline(x=p_gcn, color='b') # GCN
    plt.axvline(x=p_gan, color='g') # GAN
    plt.axvline(x=p_gin, color='k') # GIN
    plt.axvline(x=p_svr, color='m') # SVR-WL

    plt.legend(['MCQD', 'XGB', 'GCN', 'GAN', 'GIN', 'SVR-WL'])

    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()


if __name__ == "__main__":
    
    #plot_distribution(load_Tlimits(['../datasets/docking_train_tlimits.csv']), plot_name='Tlimit_dist_docking')
    #plot_distribution(load_Tlimits(['../datasets/product_train_tlimits.csv']), plot_name='Tlimit_dist_product')
    #plot_distribution(load_Tlimits(['../datasets/rand_train_tlimits.csv']), plot_name='Tlimit_dist_rand')
    

    #plot_time_dist('../datasets/dist_docking.csv', plot_name='time_dist_docking4')
    #plot_time_dist('../datasets/dist_product.csv', plot_name='time_dist_product')
    #plot_time_dist('../datasets/dist_rand.csv', plot_name='time_dist_rand')
    #plot_time_dist('../datasets/dist_dimacs.csv', plot_name='time_dist_dimacs')
    #plot_time_dist('../datasets/dist_dimacs.csv', plot_name='time_dist_MANN-a27')

    plot_time_dist_models('rand', '../datasets/rand_dist_20.csv', 'rand_dist_p', 20) # rand
    plot_time_dist_models('dense', '../datasets/dense_dist_7.csv', 'dense_dist_p', 7) # dense
    plot_time_dist_models('protein', '../datasets/protein_dist_5.csv', 'protein_dist_p', 5) # protein
    plot_time_dist_models('product', '../datasets/1_product_dist_0.csv', 'product_dist_p', 0) # product
    plot_time_dist_models('docking', '../datasets/docking_dist_7.csv', 'docking_dist_p', 7) # docking
    plot_time_dist_models('dimacs', '../datasets/dimacs_dist_0.csv', 'dimacs_dist_p', 0) # dimacs




    # train datasets: rand, product, docking
    # test datasets: rand, dense, product, docking, protein, dimacs

    '''
    datasets_test = ['rand', 'docking', 'dense']
    model_names = ['gcn', 'xgb', 'svr']
    for model_name in model_names:
        for dataset in datasets_test:
            compare_test_results(dataset, model_name)
    '''




