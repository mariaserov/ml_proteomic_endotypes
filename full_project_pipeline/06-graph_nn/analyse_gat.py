import torch
import numpy as np
import pandas as pd
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, Counter
import plotly.graph_objects as go
import plotly.express as px
from scipy.sparse import csr_matrix
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import gc

def load_attention_data(output_dir, n_folds=5):
    """Load attention weights and related data from all folds"""
    attention_dir = os.path.join(output_dir, 'attention_weights')
    attention_data = []
    
    for fold in range(1, n_folds + 1):
        try:
            data = np.load(os.path.join(attention_dir, f'attention_weights_fold{fold}.npy'),
                          allow_pickle=True).item()
            attention_data.append(data)
            print(f"Loaded attention data from fold {fold}")
        except Exception as e:
            print(f"Error loading fold {fold}: {str(e)}")
    
    return attention_data

def analyze_attention_weights(attention_data, output_dir):
    """Analyze attention weights across folds, combining bidirectional edges"""
    analysis_dir = os.path.join(output_dir, 'analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Get protein names from first fold
    protein_names = attention_data[0]['protein_names']
    
    # Aggregate attention weights across folds
    edge_attention = defaultdict(list)
    
    for data in attention_data:
        edge_index = data['edge_index']
        weights = data['attention_weights']['layer_0']['alpha']
        
        # Calculate node-wise normalisation factors
        node_degrees = defaultdict(int)
        for idx in range(edge_index.shape[1]):
            src = edge_index[0, idx]
            node_degrees[src] += 1
        
        for idx in range(0, edge_index.shape[1], 2):
            src1, dst1 = edge_index[0, idx], edge_index[1, idx]
            src2, dst2 = edge_index[0, idx+1], edge_index[1, idx+1]
            
            edge_key = tuple(sorted([
                protein_names[src1],
                protein_names[dst1]
            ]))
            
            # Raw attention weights (already softmaxed)
            weight1 = weights[idx]
            weight2 = weights[idx+1]
            
            # Normalize by source node degree
            norm_weight1 = weights[idx] * node_degrees[src1]
            norm_weight2 = weights[idx+1] * node_degrees[src2]
            
            combined_weight = (norm_weight1 + norm_weight2) / 2
            edge_attention[edge_key].append(combined_weight)
    
    # Calculate statistics
    attention_stats = {
        edge: {
            'mean': np.mean(weights),
            'std': np.std(weights),
            'source': edge[0],
            'target': edge[1]
        }
        for edge, weights in edge_attention.items()
    }
    
    # Create DataFrame of top interactions
    sorted_edges = sorted(attention_stats.items(), 
                         key=lambda x: x[1]['mean'], 
                         reverse=True)
    
    top_interactions = pd.DataFrame([
        {
            'source': stats['source'],
            'target': stats['target'],
            'mean_attention': stats['mean'],
            'std_attention': stats['std'],
            'edge_key': edge[0] + '_' + edge[1]
        }
        for edge, stats in sorted_edges[:100]
    ])
    
    # Save results
    top_interactions.to_csv(os.path.join(analysis_dir, 'top_interactions_combined.csv'), 
                           index=False)
    
    return top_interactions

def visualize_attention_weights(top_interactions, output_dir):
    """Create static visualisations of normalised attention weights"""
    viz_dir = os.path.join(output_dir, 'visualisations')
    os.makedirs(viz_dir, exist_ok=True)
    
    # 1. Network Plot
    def create_network_plot(top_n=500):
        G = nx.Graph()
        
        # Add edges with weights from averaged fold data
        for _, row in top_interactions.head(top_n).iterrows():
            G.add_edge(row['source'], row['target'], 
                      weight=float(row['mean_attention']),
                      std=float(row['std_attention']))
        
        # Calculate weighted degrees
        weighted_degree = defaultdict(float)
        for node in G.nodes():
            for neighbor in G[node]:
                weighted_degree[node] += G[node][neighbor]['weight']
        
        plt.figure(figsize=(20, 20))
        pos = nx.fruchterman_reingold_layout(G, k=5.0, iterations=200)
        
        # Draw edges with width based on weight
        max_weight = max(d['weight'] for _, _, d in G.edges(data=True))
        edge_widths = [G[u][v]['weight'] / max_weight * 3 for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=edge_widths)
        
        # Draw nodes
        max_degree = max(weighted_degree.values())
        node_sizes = [100 * weighted_degree[node] / max_degree * 5 for node in G.nodes()]
        nodes = nx.draw_networkx_nodes(G, pos,
                                     node_size=node_sizes,
                                     node_color=list(weighted_degree.values()),
                                     cmap=plt.cm.RdYlBu_r)
        
        # Calculate significance threshold using 90th percentile
        significance_threshold = np.percentile(list(weighted_degree.values()), 90)
        
        # Label significant nodes (above 90th percentile)
        significant_labels = {node: node for node in G.nodes() 
                            if weighted_degree[node] > significance_threshold}
        nx.draw_networkx_labels(G, pos, significant_labels, font_size=8,
                              bbox=dict(facecolor='white', alpha=0.7))
        
        plt.colorbar(nodes, label='Normalised Weighted Degree')
        plt.title(f'Degree-Normalised Attention Weights\n(Top {top_n} Edges, 90th percentile threshold)')
        plt.axis('off')
        plt.tight_layout()
        
        return G, weighted_degree, significance_threshold
    
    # Generate network visualisation
    print("Creating network plot...")
    G, weighted_degree, significance_threshold = create_network_plot()
    plt.savefig(os.path.join(viz_dir, 'protein_network.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Analysis Plots
    print("Creating analysis plots...")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
    
    # Attention weight distribution
    sns.histplot(data=top_interactions, x='mean_attention', bins=50, ax=ax1)
    ax1.set_title('Distribution of Normalised Attention Weights')
    ax1.set_xlabel('Normalised Attention Weight')
    ax1.set_ylabel('Count')
    
    # Attention stability (mean vs std)
    sns.scatterplot(data=top_interactions, x='mean_attention', 
                   y='std_attention', ax=ax2)
    ax2.set_title('Attention Weight Stability')
    ax2.set_xlabel('Mean Attention Weight')
    ax2.set_ylabel('Standard Deviation')
    
    # Top 20 strongest interactions
    top_20 = top_interactions.head(20).copy()
    top_20['interaction'] = top_20['source'] + ' - ' + top_20['target']
    sns.barplot(data=top_20, y='interaction', x='mean_attention', ax=ax3)
    ax3.set_title('Top 20 Strongest Interactions')
    ax3.set_xlabel('Normalised Attention Weight')
    ax3.set_ylabel('Protein Pair')
    
    plt.tight_layout()
    plt.savefig(os.path.join(viz_dir, 'attention_analysis.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print statistics about the network
    print("\nNetwork Statistics:")
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    print(f"90th percentile threshold: {significance_threshold:.4f}")
    
    print("\nTop 10 nodes by weighted degree:")
    for node, degree in sorted(weighted_degree.items(), 
                             key=lambda x: x[1], reverse=True)[:10]:
        print(f"{node}: {degree:.4f}")

def analyze_node_importance(top_interactions):
    """Analyze node importance with proper normalisation"""
    node_stats = defaultdict(lambda: {
        'incoming_sum': 0.0,
        'outgoing_sum': 0.0,
        'max_attention': 0.0,
        'edge_count': 0,
        'avg_attention': 0.0
    })
    
    # First collect all stats
    for _, row in top_interactions.iterrows():
        src, dst = row['source'], row['target']
        attention = float(row['mean_attention'])  # Ensure float
        
        # Update stats
        node_stats[src]['outgoing_sum'] += attention
        node_stats[src]['edge_count'] += 1
        node_stats[src]['max_attention'] = max(
            node_stats[src]['max_attention'], 
            attention
        )
        
        node_stats[dst]['incoming_sum'] += attention
        node_stats[dst]['edge_count'] += 1
        node_stats[dst]['max_attention'] = max(
            node_stats[dst]['max_attention'], 
            attention
        )
    
    # Then calculate normalised metrics
    for node in node_stats:
        count = max(node_stats[node]['edge_count'], 1)  # Avoid division by zero
        node_stats[node]['avg_attention'] = (
            node_stats[node]['incoming_sum'] + node_stats[node]['outgoing_sum']
        ) / (2 * count)  # Normalise by total connections
    
    # Convert to DataFrame
    node_df = pd.DataFrame.from_dict(
        node_stats, 
        orient='index'
    ).reset_index().rename(columns={'index': 'protein'})
    
    return node_df.sort_values('avg_attention', ascending=False)

def visualize_training_metrics(output_dir):
    """Plot training metrics from saved results"""
    results_dir = os.path.join(output_dir, 'results')
    viz_dir = os.path.join(output_dir, 'visualisations')
    
    try:
        with open(os.path.join(results_dir, 'cv_results.json'), 'r') as f:
            cv_results = json.load(f)['cv_results']
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
        
        for fold_idx, fold_data in enumerate(cv_results):
            epochs = range(1, len(fold_data['train_loss']) + 1)
            
            ax1.plot(epochs, fold_data['train_loss'], label=f'Fold {fold_idx+1}')
            ax2.plot(epochs, fold_data['val_auc'], label=f'Fold {fold_idx+1}')
            ax3.plot(epochs, fold_data['val_aupr'], label=f'Fold {fold_idx+1}')
        
        ax1.set_title('Training Loss')
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('Loss')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        ax2.set_title('Validation AUC')
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('AUC')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        ax3.set_title('Validation AUPR')
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('AUPR')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(viz_dir, 'training_metrics.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error plotting training metrics: {str(e)}")

if __name__ == "__main__":
    try:
        # Update base path
        base_dir = "/rds/general/project/hda_24-25/live/TDS/Group09/graph_nn"
        output_dir = os.path.join(base_dir, "gat_output_is")
        print("Starting GAT analysis...")
        
        # Load attention data
        try:
            attention_data = load_attention_data(output_dir)
        except Exception as e:
            print(f"Error loading attention data: {str(e)}")
            raise
        
        # Analyse attention weights
        print("Analysing attention weights...")
        top_interactions = analyze_attention_weights(attention_data, output_dir)
        
        # Create visualisations
        print("Creating visualisations...")
        visualize_attention_weights(top_interactions, output_dir)
        visualize_training_metrics(output_dir)
        
        # Analyse node importance
        print("Analysing node importance...")
        node_importance = analyze_node_importance(top_interactions)
        node_importance.to_csv(os.path.join(output_dir, 'analysis', 'node_importance.csv'), 
                             index=False)
        
        # After using large variables, you might want to delete them explicitly
        del attention_data
        gc.collect()
        
        print("Analysis complete!")
        
    except Exception as e:
        print(f"Error occurred during analysis: {str(e)}")
