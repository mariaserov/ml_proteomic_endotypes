import torch
import torch_geometric
import numpy as np
import pandas as pd
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from sklearn.model_selection import StratifiedKFold
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv, global_mean_pool
import os
import json
import yaml
from sklearn.metrics import roc_auc_score, average_precision_score
from torch.optim import Adam
from datetime import datetime
from collections import defaultdict

class ProteinGraphDataset:
    def __init__(self, data_path, adj_matrix_path, confidence_threshold=400):
        # Load protein expression data
        self.data = pd.read_csv(data_path)
        
        # Print initial data types and ranges for global features
        print("\nGlobal features before preprocessing:")
        print("Sex unique values:", self.data['sex.0.0'].unique())
        print("Age range:", self.data['age_at_recruitment.0.0'].min(), 
              "to", self.data['age_at_recruitment.0.0'].max())
        print("Smoking status unique values:", self.data['smoking_status_clean'].unique())
        print("Diabetes unique values:", self.data['diabetes_clean'].unique())
        
        # Verify sex is binary and convert to float
        assert set(self.data['sex.0.0'].unique()).issubset({0, 1}), \
            "Sex values should only be 0 or 1"
        self.data['sex.0.0'] = self.data['sex.0.0'].astype(float)
        
        # Normalize age
        age_mean = self.data['age_at_recruitment.0.0'].mean()
        age_std = self.data['age_at_recruitment.0.0'].std()
        self.data['age_at_recruitment.0.0'] = (self.data['age_at_recruitment.0.0'] - age_mean) / age_std
        
        # One-hot encode smoking status (3 levels: 0, 1, 2)
        smoking_dummies = pd.get_dummies(self.data['smoking_status_clean'], prefix='smoking')
        smoking_dummies = smoking_dummies.astype(float)  # Convert boolean to float
        self.data = pd.concat([self.data, smoking_dummies], axis=1)
        
        # Convert diabetes to float
        self.data['diabetes_clean'] = self.data['diabetes_clean'].astype(float)
        
        # Store normalization parameters
        self.global_feature_params = {
            'age_mean': age_mean,
            'age_std': age_std
        }
        
        # Get protein columns
        self.protein_cols = [col for col in self.data.columns if col not in [
            'Row.names', 'sex.0.0', 'age_at_recruitment.0.0', 'smoking_status_clean',
            'smoking_0', 'smoking_1', 'smoking_2', 'diabetes_clean', 'incident_case'
        ]]
        
        # Verify protein count
        assert len(self.protein_cols) == 1456
        
        # Normalize protein values
        self.data[self.protein_cols] = self.data[self.protein_cols].astype(float)
        self.protein_means = self.data[self.protein_cols].mean()
        self.protein_stds = self.data[self.protein_cols].std()
        self.data[self.protein_cols] = (self.data[self.protein_cols] - self.protein_means) / self.protein_stds
        
        # Save normalization parameters
        normalization_params = {
            'means': self.protein_means.to_dict(),
            'stds': self.protein_stds.to_dict()
        }
        with open(os.path.join(os.path.dirname(data_path), 'normalization_params.json'), 'w') as f:
            json.dump(normalization_params, f)
        
        # Process adjacency matrix
        self.adj_matrix = np.load(adj_matrix_path)
        assert self.adj_matrix.shape == (1456, 1456)
        
        # Create edge list
        rows, cols = np.where(self.adj_matrix >= confidence_threshold)
        edge_list = []
        edge_weights = []
        
        for r, c in zip(rows, cols):
            if r != c:
                edge_list.extend([[r, c], [c, r]])
                edge_weights.extend([self.adj_matrix[r, c], self.adj_matrix[r, c]])
        
        # Convert to PyTorch tensors
        self.edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        self.edge_weights = torch.tensor(edge_weights, dtype=torch.float)
        
        # Normalize edge weights
        self.edge_weights = (self.edge_weights - self.edge_weights.min()) / \
                          (self.edge_weights.max() - self.edge_weights.min())
        
        # Calculate pos_weight for loss function
        num_positives = self.data['incident_case'].sum()
        num_negatives = len(self.data) - num_positives
        self.pos_weight = num_negatives / num_positives
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        features = torch.tensor(self.data.iloc[idx][self.protein_cols].values, 
                              dtype=torch.float).unsqueeze(1)
        
        # Get global features (sex, age, smoking status (3 levels), diabetes)
        global_features = torch.tensor(
            self.data.iloc[idx][['sex.0.0', 'age_at_recruitment.0.0', 
                                'smoking_0', 'smoking_1', 'smoking_2',
                                'diabetes_clean']].values,
            dtype=torch.float
        ).unsqueeze(0)
        
        label = torch.tensor(self.data.iloc[idx]['incident_case'], 
                           dtype=torch.float)
        
        graph = Data(
            x=features,
            edge_index=self.edge_index.clone(),
            edge_attr=self.edge_weights.clone(),
            global_features=global_features,
            y=label
        )
        
        graph.num_nodes = features.size(0)
        return graph

class ProteinGAT(nn.Module):
    def __init__(self, num_proteins, protein_features, global_features,
                 hidden_dim=64, num_layers=2, dropout=0.1, edge_dropout=0.05):
        super(ProteinGAT, self).__init__()
        
        self.num_layers = num_layers
        self.dropout = dropout
        self.edge_dropout = edge_dropout
        
        self.protein_embedding = nn.Linear(protein_features, hidden_dim)
        self.global_embedding = nn.Linear(6, hidden_dim)  # Updated for 6 global features
        
        self.gat_layers = nn.ModuleList()
        self.gat_layers.append(
            GATConv(hidden_dim, hidden_dim, dropout=dropout,
                   add_self_loops=False, edge_dim=1)
        )
        
        for _ in range(num_layers - 1):
            self.gat_layers.append(
                GATConv(hidden_dim, hidden_dim, dropout=dropout,
                       add_self_loops=False, edge_dim=1)
            )
        
        self.combine_layer = nn.Linear(hidden_dim * 2, hidden_dim)
        self.output_layer = nn.Linear(hidden_dim, 1)
        
        # Store attention weights
        self.attention_weights = {}
        
    def forward(self, data):
        x, edge_index, edge_attr, global_features, batch = (
            data.x,
            data.edge_index,
            data.edge_attr.unsqueeze(1),
            data.global_features,
            data.batch
        )
        
        x = self.protein_embedding(x)
        global_feat = self.global_embedding(global_features)
        
        if self.training and self.edge_dropout > 0:
            edge_mask = torch.rand(edge_index.size(1), device=edge_index.device) > self.edge_dropout
            edge_index = edge_index[:, edge_mask]
            edge_attr = edge_attr[edge_mask]
        
        self.attention_weights.clear()
        for i, gat_layer in enumerate(self.gat_layers):
            x, (edge_index, alpha) = gat_layer(x, edge_index, edge_attr=edge_attr, 
                                             return_attention_weights=True)
            
            if not self.training:
                self.attention_weights[f'layer_{i}'] = {
                    'edge_index': edge_index,
                    'alpha': alpha,
                    'edge_attr': edge_attr
                }
            
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        graph_repr = global_mean_pool(x, batch)
        combined = torch.cat([graph_repr, global_feat], dim=1)
        combined = self.combine_layer(combined)
        combined = F.relu(combined)
        
        return self.output_layer(combined)

def train_epoch(model, loader, optimizer, device, pos_weight):
    model.train()
    total_loss = 0
    
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        
        logits = model(batch).squeeze(1)
        loss = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight]).to(device))(
            logits, batch.y.float())
        
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item() * batch.num_graphs
    
    return total_loss / len(loader.dataset)

def evaluate(model, loader, device):
    model.eval()
    predictions = []
    targets = []
    all_attention_weights = defaultdict(lambda: defaultdict(list))
    
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            predictions.append(torch.sigmoid(out).cpu())
            targets.append(batch.y.cpu())
            
            # Store attention weights from this batch
            for layer_name, layer_weights in model.attention_weights.items():
                all_attention_weights[layer_name]['edge_index'].append(layer_weights['edge_index'].cpu())
                all_attention_weights[layer_name]['alpha'].append(layer_weights['alpha'].cpu())
                all_attention_weights[layer_name]['edge_attr'].append(layer_weights['edge_attr'].cpu())
            
            model.attention_weights.clear()  # Clear for next batch
    
    # Combine attention weights across all batches
    combined_attention = {}
    for layer_name, layer_data in all_attention_weights.items():
        combined_attention[layer_name] = {
            'edge_index': torch.cat(layer_data['edge_index'], dim=1),
            'alpha': torch.cat(layer_data['alpha'], dim=0),
            'edge_attr': torch.cat(layer_data['edge_attr'], dim=0)
        }
    
    # Store combined attention weights
    model.attention_weights = combined_attention
    
    predictions = torch.cat(predictions, dim=0)
    targets = torch.cat(targets, dim=0)
    
    return predictions, targets

def train_model(dataset, output_dir, n_splits=5, epochs=200, batch_size=32,
                learning_rate=0.0005, hidden_dim=64, num_layers=2,
                dropout=0, edge_dropout=0, patience=15,
                device='cuda' if torch.cuda.is_available() else 'cpu'):
    
    # Create output directories
    model_dir = os.path.join(output_dir, 'models')
    results_dir = os.path.join(output_dir, 'results')
    attention_dir = os.path.join(output_dir, 'attention_weights')
    
    for dir_path in [model_dir, results_dir, attention_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Save configuration
    config = {
        'n_splits': n_splits,
        'epochs': epochs,
        'batch_size': batch_size,
        'learning_rate': learning_rate,
        'hidden_dim': hidden_dim,
        'num_layers': num_layers,
        'dropout': dropout,
        'edge_dropout': edge_dropout,
        'patience': patience,
        'device': device,
        'timestamp': datetime.now().strftime('%Y%m%d_%H%M%S')
    }
    
    with open(os.path.join(output_dir, 'config.yaml'), 'w') as f:
        yaml.dump(config, f)
    
    # Initialize cross-validation
    kfold = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    cv_results = []
    
    all_labels = [data.y.item() for data in dataset]
    
    for fold, (train_idx, val_idx) in enumerate(kfold.split(range(len(dataset)), all_labels)):
        print(f"\nTraining Fold {fold + 1}/{n_splits}")
        
        train_loader = DataLoader([dataset[i] for i in train_idx], batch_size=batch_size, shuffle=True)
        val_loader = DataLoader([dataset[i] for i in val_idx], batch_size=batch_size)
        
        model = ProteinGAT(
            num_proteins=1456,
            protein_features=1,
            global_features=6,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            dropout=dropout,
            edge_dropout=edge_dropout
        ).to(device)
        
        optimizer = Adam(model.parameters(), lr=learning_rate)
        
        best_val_auc = 0
        best_epoch = 0
        no_improve = 0
        fold_results = {'train_loss': [], 'val_auc': [], 'val_aupr': []}
        
        for epoch in range(epochs):
            train_loss = train_epoch(model, train_loader, optimizer, device, dataset.pos_weight)
            
            # Evaluate and get attention weights for full validation set
            val_pred, val_true = evaluate(model, val_loader, device)
            val_auc = roc_auc_score(val_true, val_pred)
            val_aupr = average_precision_score(val_true, val_pred)
            
            fold_results['train_loss'].append(train_loss)
            fold_results['val_auc'].append(val_auc)
            fold_results['val_aupr'].append(val_aupr)
            
            print(f"Epoch {epoch+1}: Loss = {train_loss:.4f}, Val AUC = {val_auc:.4f}, Val AUPR = {val_aupr:.4f}")
            
            if val_auc > best_val_auc:
                best_val_auc = val_auc
                best_epoch = epoch
                no_improve = 0
                
                # Save model and attention weights (now containing full validation set)
                torch.save(model.state_dict(), 
                         os.path.join(model_dir, f'best_model_fold{fold+1}.pt'))
                
                attention_data = {
                    'edge_index': dataset.edge_index.cpu().numpy(),
                    'attention_weights': {
                        k: {
                            'edge_index': v['edge_index'].cpu().numpy(),
                            'alpha': v['alpha'].cpu().numpy(),
                            'edge_attr': v['edge_attr'].cpu().numpy()
                        }
                        for k, v in model.attention_weights.items()
                    },
                    'protein_names': dataset.protein_cols,
                    'validation_indices': val_idx  # Add this to track which samples were used
                }
                np.save(os.path.join(attention_dir, f'attention_weights_fold{fold+1}.npy'),
                        attention_data)
            else:
                no_improve += 1
                if no_improve >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break
        
        fold_results['best_epoch'] = best_epoch
        fold_results['best_val_auc'] = best_val_auc
        cv_results.append(fold_results)
        
        # Clear memory
        del model
        torch.cuda.empty_cache()
    
    # Save overall results
    overall_results = {
        'mean_best_val_auc': np.mean([r['best_val_auc'] for r in cv_results]),
        'std_best_val_auc': np.std([r['best_val_auc'] for r in cv_results]),
        'cv_results': cv_results
    }
    
    with open(os.path.join(results_dir, 'cv_results.json'), 'w') as f:
        json.dump(overall_results, f, indent=4)
    
    return overall_results

if __name__ == "__main__":
    try:
        # Update base path
        base_dir = "/rds/general/project/hda_24-25/live/TDS/Group09/graph_nn"
        output_dir = os.path.join(base_dir, "gat_output_is")
        print("Starting GAT training...")
        
        dataset = ProteinGraphDataset(
            data_path=os.path.join(base_dir, "data/gat_input_data_is.csv"),
            adj_matrix_path=os.path.join(base_dir, "output/protein_adjacency_matrix.npy"),
            confidence_threshold=400
        )
        
        results = train_model(dataset, output_dir)
        print(f"Training complete! Mean CV AUC: {results['mean_best_val_auc']:.4f} ± {results['std_best_val_auc']:.4f}")
        
    except Exception as e:
        print(f"Error occurred during execution: {str(e)}")
