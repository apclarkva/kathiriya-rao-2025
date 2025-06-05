import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import numpy as np



def visualize_sub(pd_net, tf_genes, targets, name):
    
    string = ''
    for s in targets:
        string += f'{s}, '

    targets += tf_genes

    net, G1, G2 = plot_subnetwork(pd_net, tf_genes, targets)

    net.show_buttons(filter_=['physics'])
    net.show(f'data/{name}', notebook=False)


def plot_subnetwork(network, tf_genes, targets):
    graph_df = {'source': [],
                 'target': [],
                 'weight': [],
                 'edge_cols': []}

    for curr_gene in tf_genes:
        curr_gene_df = network[network['source'] == curr_gene]
        gene_df_targs = curr_gene_df[curr_gene_df['target'].isin(targets)]

        for index, row in gene_df_targs.iterrows():
            graph_df['source'].append(row['source'])
            graph_df['target'].append(row['target'])
            graph_df['weight'].append(row['coef_mean'])
            if row['coef_mean'] > 0:
                graph_df['edge_cols'].append('green')
            else:
                graph_df['edge_cols'].append('red')

    graph_df = pd.DataFrame(graph_df)

    G = nx.from_pandas_edgelist(graph_df, source='source',
                                    target='target',
                                    edge_attr='weight',
                                    create_using=nx.DiGraph())

    from pyvis.network import Network
    net = Network(notebook=True, directed=True, height='500px', width='1000px')
    #net.force_atlas_2based(gravity=-10, spring_length=80)

    #For 12,000 edge networks
    #net.force_atlas_2based(gravity=-10, spring_length=40)

    #For priority genes
    net.force_atlas_2based(gravity=-30, spring_strength=.005)


    net.from_nx(G)

    for edge in net.edges:
        if edge['width'] > 0:
            edge['color'] = 'grey'
            edge['arrowStrikethrough'] = False
        else:
            edge['color'] = 'red'
            edge['arrows'] = {"to": {"enabled": True, "type": "bar",
                                     "scaleFactor":1.4}}

        edge['value'] = np.abs(edge['width'])


    for node in net.nodes:
        if node['id'] in tf_genes:
            #node['shape'] = 'circle'
            node['color'] = 'darkgrey'
            node['size'] = 20 
            node['font'] = {'size': 50,
                            'color': 'black',
                            'style': 'italic',
                            'multi': 'html',
                            'vadjust': -30 # Move text closer to the node (negative values pull it up)
            }
            # Settings for 12,000 edge
            ##node['shape'] = 'circle'
            #node['color'] = 'darkgrey'
            #node['size'] = 40 
            #node['font'] = {'size': 125,
            #                'color': 'black',
            #                'style': 'italic',
            #                'multi': 'html',
            #                'vadjust': -70 # Move text closer to the node (negative values pull it up)
            #}
        else:
            node['color'] = 'darkgrey'
            node['font'] = {'size': 2,
                            'color': 'black',
                            'style':'italic',
                            'multi': 'html'}
            node['size'] = 6 

        if node['id'] in ['Gata4', 'Nr2f2']:
            node['shape'] = 'diamond'
            node['color'] = 'darkgrey'
            node['size'] = 20
            node['font'] = {'size': 50,
                            'color': 'black',
                            'multi': 'html'}

        node['label'] = f'<i>{node["id"]}</i>'



    G_act = nx.from_pandas_edgelist(graph_df[graph_df.weight > 0], source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    G_inact = nx.from_pandas_edgelist(graph_df[graph_df.weight < 0], source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())


    return net, G_act, G_inact


def get_targets(net, genes):
    c = []
    c = [c + list(net[net.source == curr_tf].target.values) for curr_tf in genes]
    targets = [item for sublist in c for item in sublist]
    return targets


def main():
    #atr_net = pd.read_csv('./data/atrial-vent-040125/atrial_12000.csv')
    #vent_net = pd.read_csv('./data/atrial-vent-040125/vent_12000.csv')

    #genes_for_subnetwork = ['TBX5', 'NR2F2', 'MEF2C', 'NKX2-5', 'MEIS2', 'GATA4', 'HAND1', 'IRX4', 'HEY2', 'PLAGL1', 'NR2F1', 'HES4', 'E2F1', 'HEY1']
    
    #visualize_sub(atr_net, genes_for_subnetwork, get_targets(atr_net, genes_for_subnetwork), 'atr_12000.html')
    #visualize_sub(vent_net, genes_for_subnetwork, get_targets(vent_net, genes_for_subnetwork), 'vent_12000.html')

    #genes_for_subnetwork = ['TBX5', 'NR2F2', 'MEF2C', 'NKX2-5', 'MEIS2', 'GATA4', 'HAND1', 'IRX4', 'HEY2', 'PLAGL1', 'NR2F1', 'HES4', 'E2F1', 'HEY1', 'TCF21', 'TEAD1', 'TCF12']

    wt_net = pd.read_csv('./data/genotype_2500402/wt_12000.csv')
    het_net = pd.read_csv('./data/genotype_2500402/het_12000.csv')
    ko_net = pd.read_csv('./data/genotype_2500402/ko_12000.csv')

    #visualize_sub(wt_net, genes_for_subnetwork, get_targets(wt_net, genes_for_subnetwork), 'wt_12000.html')
    #visualize_sub(het_net, genes_for_subnetwork, get_targets(het_net, genes_for_subnetwork), 'het_12000.html')
    #visualize_sub(ko_net, genes_for_subnetwork, get_targets(ko_net, genes_for_subnetwork), 'ko_12000.html')

    priority_genes = pd.read_excel("data/genes_for_subnetwork-250316.xlsx", sheet_name=0).values
    priority_genes = priority_genes[~pd.isna(priority_genes)]

    priority_genes = list(priority_genes) + ['TBX5'] 
    visualize_sub(wt_net, priority_genes, priority_genes, 'chd_af_wt.html')
    visualize_sub(het_net, priority_genes, priority_genes, 'chd_af_het.html')
    visualize_sub(ko_net, priority_genes, priority_genes, 'chd_af_ko.html')

    


if __name__ == '__main__':
    main()
