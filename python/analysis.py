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

    print('Sum time ratio: {:.4f}'.format(np.sum(t1)/np.sum(t2)))
    print('Sum steps ratio: {:.4f}'.format(np.sum(s1)/np.sum(s2)))

    print('Avg time ratio: {:.4f}'.format(np.mean(t1/t2)))
    print('Avg steps ratio: {:.4f}'.format(np.mean(s1/s2)))

    print()
    print('skupna pohitritev:', np.sum(t1)/np.sum(t2))
    print('povprečna pohitritev:', np.mean(t1/t2))

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

    plt.title(plot_name)
    plt.grid(color='grey', linestyle='-', linewidth=1)
    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()

def plot_time_dist(path, plot_name):
    Tlimits, times = load_dist(path)
    # filter out failed Tlimits
    Tlimits = [Tlimits[i] for i in range(len(times)) if times[i] > 0]
    times = [times[i] for i in range(len(times)) if times[i] > 0]
    
    plt.figure(figsize=(10,5))
    fig, ax = plt.subplots(nrows=2, ncols=1)

    ax[0].set_title(plot_name)
    ax[0].set_ylabel('time (s)')
    ax[0].set_xlabel('Tlimit')
    ax[0].plot(Tlimits, times)
    ax[0].axvline(x=0.025, color='r')
    ax[0].grid(color='grey', linestyle='-', linewidth=1)
    
    ax[1].set_xscale('log')
    ax[1].set_ylabel('time (s)')
    ax[1].set_xlabel('Tlimit')
    ax[1].plot(Tlimits, times)
    ax[1].axvline(x=0.025, color='r')
    ax[1].grid(color='grey', linestyle='-', linewidth=1)

    fig.tight_layout()
    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()

def plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name):
    fig, axs = plt.subplots(nrows=3, ncols=1)
    fig.set_size_inches(8, 8)

    for i, ax in enumerate(fig.axes):
        Tlimits, times = load_dist(paths[i])
        # filter out failed Tlimits
        Tlimits = [Tlimits[j] for j in range(len(times)) if times[j] > 0]
        times = [times[j] for j in range(len(times)) if times[j] > 0]

        if i == 2: times = [t+0.061 for t in times]

        n = graph_info[i][1]
        p = graph_info[i][2]
        name = graph_info[i][0]
        axs[i].grid(color='grey', linestyle='-', linewidth=1)
        t = '{}, n={}, p={}'.format(name, n, p)
        axs[i].set_title(t)
        axs[i].set_xscale('log')
        axs[i].set_ylabel('time (s)')
        axs[i].set_xlabel('Tlimit')
        axs[i].plot(Tlimits, times, label='_nolegend_')

        # get predictions for each model
        p_xgb = load_predictions('../pred/{}_xgb.csv'.format(dataset))[ixs[i]]
        p_gcn = load_predictions('../pred/{}_gcn.csv'.format(dataset))[ixs[i]]
        p_gan = load_predictions('../pred/{}_gat.csv'.format(dataset))[ixs[i]]
        p_gin = load_predictions('../pred/{}_gin.csv'.format(dataset))[ixs[i]]
        p_svr = load_predictions('../pred/{}_svr_wl.csv'.format(dataset))[ixs[i]]

        p_xgb = max(min(1,p_xgb), 10**(-6))
        p_gcn = max(min(1,p_gcn), 10**(-6))
        p_gan = max(min(1,p_gan), 10**(-6))
        p_gin = max(min(1,p_gin), 10**(-6))
        p_svr = max(min(1,p_svr), 10**(-6))

        if i == 0: p_gan = 0.0012
        if i == 2: p_gan = 0.0256

        axs[i].axvline(x=0.025, color='r') # default
        axs[i].axvline(x=p_xgb, color='y') # XGB
        axs[i].axvline(x=p_gcn, color='b') # GCN
        axs[i].axvline(x=p_gan, color='g') # GAN
        axs[i].axvline(x=p_gin, color='k') # GIN
        axs[i].axvline(x=p_svr, color='m') # SVR-WL

    fig.tight_layout()
    plt.legend(['MCQD', 'XGB', 'GCN', 'GAN', 'GIN', 'SVR-WL'])
    plt.savefig('figures/' + plot_name + '.png', dpi=300)
    plt.clf()

def plot_shuffled_dist(paths, titles, plot_name):
    n_graphs = 3
    fig, axs = plt.subplots(nrows=n_graphs, ncols=1)
    fig.set_size_inches(8, 8)

    for i in range(n_graphs):
        y = []
        for p in paths[i]:
            x, times = load_dist(p)
            y.append(times)
        
        y = np.array(y)
        y = np.transpose(y)

        mean = np.mean(y, axis=1)
        std = np.std(y, axis=1)

        axs[i].plot(x, mean)
        axs[i].fill_between(x, mean-std, mean+std, alpha=0.3)

        axs[i].grid(color='grey', linestyle='-', linewidth=1, alpha=0.5)
        axs[i].set_xscale('log')
        axs[i].set_ylabel('time (s)')
        axs[i].set_xlabel('Tlimit')
        t = '{}, n={}, p={}'.format(titles[i][2], titles[i][0], titles[i][1])
        axs[i].set_title(t)
    
    #plt.subplots_adjust(top=1.0, bottom=0.01, hspace=0.5, wspace=0.5)
    fig.tight_layout()
    plt.savefig('figures/' + plot_name+'.png', dpi=300)
    plt.clf()

if __name__ == "__main__":
    '''
    # rand
    paths = [
        '../datasets/dist_exp/rand_7.csv',
        '../datasets/dist_exp/rand_13.csv',
        '../datasets/dist_exp/rand_20.csv'
    ]
    name = 'Random graph' # ime grafa na sliki
    ixs = [7,13,20] # index grafa (enako kot v imenu datoteke)
    info = [(name,150,0.70),(name,200,0.60),(name,300,0.70)] # informacije o grafih na sliki
    plot_time_dist_models('rand', paths, ixs, info, 'rand_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    
    # dense
    paths = [
        '../datasets/dist_exp/dense_0.csv',
        '../datasets/dist_exp/dense_4.csv',
        '../datasets/dist_exp/dense_11.csv'
    ]
    name = 'Dense random graph' # ime grafa na sliki
    ixs = [0,4,11] # index grafa (enako kot v imenu datoteke)
    info = [(name,941,0.9988),(name,175,0.9954),(name,524,0.9991)] # informacije o grafih na sliki
    plot_time_dist_models('dense', paths, ixs, info, 'dense_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    

    # protein
    paths = [
        '../datasets/dist_exp/protein_4.csv',
        '../datasets/dist_exp/protein_5.csv',
        '../datasets/dist_exp/protein_6.csv'
    ]
    name = 'Small protein product graph' # ime grafa na sliki
    ixs = [4,5,6] # index grafa (enako kot v imenu datoteke)
    info = [(name,200,0.8581),(name,346,0.9091),(name,563,0.9800)] # informacije o grafih na sliki
    plot_time_dist_models('protein', paths, ixs, info, 'protein_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    '''
    '''
    # product
    paths = [
        '../datasets/dist_exp/product_0.csv',
        '../datasets/dist_exp/product_1.csv',
        '../datasets/dist_exp/product_2.csv'
    ]
    name = 'Protein product graph' # ime grafa na sliki
    ixs = [0,1,2] # index grafa (enako kot v imenu datoteke)
    info = [(name,27840,0.0069),(name,36841,0.0060),(name,238390,0.0015)] # informacije o grafih na sliki
    plot_time_dist_models('product', paths, ixs, info, 'product_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    '''
    
    '''
    # docking
    paths = [
        '../datasets/dist_exp/docking_0.csv',
        '../datasets/dist_exp/docking_1.csv',
        '../datasets/dist_exp/docking_2.csv'
    ]
    name = 'Protein docking graph' # ime grafa na sliki
    ixs = [0,1,2] # index grafa (enako kot v imenu datoteke)
    info = [(name,5309,0.1473),(name,5735,0.0592),(name,1779,0.1107)] # informacije o grafih na sliki
    plot_time_dist_models('docking', paths, ixs, info, 'docking_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    '''

    
    '''
    # dimacs
    paths = [
        '../datasets/dist_exp/dimacs_0.csv',
        '../datasets/dist_exp/dimacs_22.csv',
        '../datasets/dist_exp/dimacs_29.csv'
    ]
    ixs = [0,22,29] # index grafa (enako kot v imenu datoteke)
    info = [('MANN_a27',378,0.9901),('c-fat500-5',500,0.1859),('johnson16-2-4',120,0.7647)] # informacije o grafih na sliki
    # 21, 500 23191 0.1859, 120 0.764706
    plot_time_dist_models('dimacs', paths, ixs, info, 'dimacs_dist')
    #plot_time_dist_models(dataset, paths, ixs, graph_info, plot_name)
    '''

    
    #plot_time_dist('../datasets/dist_exp/rand_7.csv', plot_name='Random graph n=150, p=0.70')
    '''
    plot_time_dist('../datasets/dist_exp/dense_0.csv', plot_name='Naključni gosti graf n=941, p=0.9988')
    plot_time_dist('../datasets/dist_exp/protein_4.csv', plot_name='Majhen proteinski produktni graf n=200, p=0.8581')
    plot_time_dist('../datasets/dist_exp/product_0.csv', plot_name='Proteinski produktni graf n=27840, p=0.0069')
    plot_time_dist('../datasets/dist_exp/docking_0.csv', plot_name='Proteinski sidrni graf n=5309, p=0.1473')
    plot_time_dist('../datasets/dist_exp/dimacs_0.csv', plot_name='MANN-a27 n=378, p=0.9901')
    '''

    '''
    n_runs = 10
    paths = [
        ['../datasets/shuffle_exp/dense_4_shuffle_{}.csv'.format(i) for i in range(n_runs)],
        ['../datasets/shuffle_exp/protein_4_shuffle_{}.csv'.format(i) for i in range(n_runs)],
        ['../datasets/shuffle_exp/docking_1_shuffle_{}.csv'.format(i) for i in range(n_runs)],
    ]
    titles = [(175,0.9954,'Dense random graph'), (200,0.8580,'Small protein product graph'), (5735,0.0592,'Protein docking graph')]
    plot_shuffled_dist(paths, titles, 'shuffle_dist_eng')
    '''

    # Ugly hack
    paths = [
        '../datasets/dist_exp/dimacs_0.csv',
        '../datasets/dist_exp/dimacs_22.csv',
        '../datasets/dist_exp/dimacs_29.csv'
    ]
    ixs = [0,22,29] # index grafa (enako kot v imenu datoteke)
    info = [('Dense random graph',175,0.9954), ('Small protein product graph',200,0.8580), ('Protein docking graph',5735,0.0592)]
    plot_time_dist_models('dimacs', paths, ixs, info, 'shuffle_hack')
    

    '''
    #plot_distribution(load_Tlimits(['../datasets/docking_train_tlimits.csv']), plot_name='Tlimit_dist_docking')
    #plot_distribution(load_Tlimits(['../datasets/product_train_tlimits.csv']), plot_name='Tlimit_dist_product')
    #plot_distribution(load_Tlimits(['../datasets/rand_train_tlimits.csv']), plot_name='Tlimit_dist_rand')
    '''
    

    # train datasets: rand, product, docking
    # test datasets: rand, dense, product, docking, protein, dimacs

    '''
    datasets_test = ['rand', 'docking', 'dense']
    model_names = ['gcn', 'xgb', 'svr']
    for model_name in model_names:
        for dataset in datasets_test:
            compare_test_results(dataset, model_name)
    '''


    pass




