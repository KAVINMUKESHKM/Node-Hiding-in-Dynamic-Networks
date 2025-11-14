import os
import sys
import networkx as nx
import numpy as np
from collections import defaultdict
import csv
import time
import argparse
import json

# Import the angel library
try:
    import angel as a
except ImportError:
    print("ERROR: Could not find the 'angel' library.")
    sys.exit(1)

# Import optional libraries for plotting
try:
    import pandas as pd
    import matplotlib.pyplot as plt
    PLOTLIBS_AVAILABLE = True
except ImportError:
    print("WARNING: 'pandas' or 'matplotlib' not found. Plotting functions will be disabled.")
    PLOTLIBS_AVAILABLE = False


class BIHDynamicHiding:
    """
    Based Importance Hiding (BIH) algorithm for dynamic networks
    with ArchAngel community detection
    """
    
    def __init__(self, base_path, threshold=0.35, match_threshold=0.35, min_comsize=3,
                 simulate_targets=True, debug_sim=False, max_candidate_sims=None,
                 enable_hillclimb=True, hillclimb_max_iter=3):
        self.base_path = base_path
        self.threshold = threshold
        self.match_threshold = match_threshold
        self.min_comsize = min_comsize
        # Controls whether to simulate target communities and pick the best
        self.simulate_targets = simulate_targets
        # When True, write per-snapshot simulation debug files into comparison dir
        self.debug_sim = debug_sim
        # Limit how many candidate communities to simulate (top-K); None = all
        self.max_candidate_sims = max_candidate_sims
        # Hill-climb options
        self.enable_hillclimb = enable_hillclimb
        self.hillclimb_max_iter = hillclimb_max_iter
        
        # Create main output directory
        self.output_dir = os.path.join(base_path, "BIH_Results")
        self.networks_dir = os.path.join(self.output_dir, "networks")
        self.comparison_dir = os.path.join(self.output_dir, "comparison")
        
        for path in [self.output_dir, self.networks_dir, self.comparison_dir]:
            os.makedirs(path, exist_ok=True)
    
    def load_snapshot_network(self, ncol_file, snapshot_id):
        """Load a specific snapshot from the ncol file"""
        G = nx.Graph()
        with open(ncol_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    node1, node2, snap = parts[0], parts[1], parts[2]
                    if snap == str(snapshot_id):
                        G.add_edge(node1, node2)
        return G
    
    def load_communities(self, comm_file):
        """Load communities from ArchAngel output file"""
        communities = {}
        with open(comm_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    comm_id = int(parts[0])
                    nodes = eval(parts[1])  # Convert string list to actual list
                    communities[comm_id] = set(nodes)
        return communities
    
    def find_overlapping_nodes(self, communities):
        """Find all nodes that belong to multiple communities"""
        node_memberships = defaultdict(list)
        for comm_id, nodes in communities.items():
            for node in nodes:
                node_memberships[node].append(comm_id)
        
        overlapping = {node: comms for node, comms in node_memberships.items() 
                       if len(comms) > 1}
        return overlapping
    
    def calculate_importance_degree(self, G, node, community_nodes):
        """
        Calculate importance degree of a node in a community
        Enhanced version that considers community density and relative connectivity
        """
        if node not in G:
            return 0
        
        degree = G.degree(node)
        if degree == 0:
            return 0
        
        neighbors_in_comm = set(G.neighbors(node)) & community_nodes
        
        # Calculate community density
        comm_subgraph = G.subgraph(community_nodes)
        possible_edges = len(community_nodes) * (len(community_nodes) - 1) / 2
        actual_edges = comm_subgraph.number_of_edges()
        density = actual_edges / possible_edges if possible_edges > 0 else 0
        
        # Calculate relative connectivity (what fraction of node's edges are in this community)
        connectivity_ratio = len(neighbors_in_comm) / degree
        
        importance_sum = 0
        for neighbor in neighbors_in_comm:
            neighbor_set = set(G.neighbors(neighbor)) & community_nodes
            common_neighbors = len(neighbors_in_comm & neighbor_set)
            importance_sum += common_neighbors
        
        # Adjust importance based on density and connectivity ratio
        raw_importance = (importance_sum * (degree - 1)) / degree if degree > 0 else 0
        
        # Scale importance:
        # - Higher for target community (want high-density regions)
        # - Lower for other communities (easier to remove from high-density regions)
        return raw_importance * (1 + density) * connectivity_ratio
    
    def get_node_with_highest_importance(self, G, target_node, candidate_nodes, community_nodes):
        """Find node with highest importance for BIH operation"""
        max_importance = -1
        best_node = None
        
        for candidate in candidate_nodes:
            if candidate == target_node:
                continue
            importance = self.calculate_importance_degree(G, candidate, community_nodes)
            if importance > max_importance:
                max_importance = importance
                best_node = candidate
        
        return best_node
    
    def try_hide_in_community(self, G, target_node, comm_id, communities, comm_metrics):
        """
        Attempt to hide a node in a specific community and return success score
        Lower score = better hiding
        """
        G_temp = G.copy()
        target_community = communities[comm_id]
        neighbors = set(G_temp.neighbors(target_node))
        
        # Calculate overlap with target community
        overlap = len(neighbors & target_community) / (len(neighbors) + 1e-6)
        
        # Try strategic rewiring
        nonneighbors = target_community - neighbors - {target_node}
        added = 0
        target_avg_degree = comm_metrics.get(comm_id, {'avg_degree': 5})['avg_degree']
        for n in sorted(nonneighbors, 
                        key=lambda x: len(set(G_temp.neighbors(x)) & target_community),
                        reverse=True)[:3]:  # Try top 3 candidates
            if added < target_avg_degree:
                G_temp.add_edge(target_node, n)
                added += 1
        
        # Calculate final score (lower is better)
        new_neighbors = set(G_temp.neighbors(target_node))
        new_overlap = len(new_neighbors & target_community) / (len(new_neighbors) + 1e-6)
        
        # Score combines overlap improvement and degree preservation
        score = (1 - new_overlap) + abs(len(new_neighbors) - len(neighbors)) / len(neighbors)
        return score, G_temp

    def hide_from_single_community(self, G, target_node, from_comm, to_comm, communities, comm_metrics):
        """
        Progressive strategy: try to hide a node from one specific community
        while strengthening its membership in the target community
        """
        G_temp = G.copy()
        from_nodes = communities[from_comm]
        to_nodes = communities[to_comm]
        
        # Phase 1: Weaken connections to source community
        neighbors = set(G_temp.neighbors(target_node))
        from_neighbors = neighbors & from_nodes
        
        # Keep minimal connections based on community density
        min_edges = max(1, int(comm_metrics[from_comm]['density'] * len(from_neighbors) * 0.2))
        
        # Sort neighbors by their importance to the from_community
        to_remove = sorted(
            from_neighbors,
            key=lambda n: (
                len(set(G_temp.neighbors(n)) & from_nodes) / G_temp.degree(n),
                -len(set(G_temp.neighbors(n)) & to_nodes)
            ),
            reverse=True
        )
        
        # Remove edges while preserving minimal connectivity
        ops = []
        for n in to_remove[:-min_edges]:
            G_temp.remove_edge(target_node, n)
            ops.append(f"Remove: {target_node}-{n}")
            
        # Phase 2: Strengthen connections to target community
        non_neighbors = to_nodes - set(G_temp.neighbors(target_node)) - {target_node}
        target_avg_degree = comm_metrics[to_comm]['avg_degree']
        
        # Sort potential new neighbors by their embeddedness in target community
        candidates = sorted(
            non_neighbors,
            key=lambda n: (
                len(set(G_temp.neighbors(n)) & to_nodes) / G_temp.degree(n),
                -len(set(G_temp.neighbors(n)) & from_nodes)
            ),
            reverse=True
        )
        
        # Add strategic connections
        added = 0
        for n in candidates:
            if added < target_avg_degree * 0.7:  # Add up to 70% of average degree
                G_temp.add_edge(target_node, n)
                ops.append(f"Add: {target_node}-{n}")
                added += 1
                
        return G_temp, ops

    def apply_bih_to_snapshot(self, G, communities, target_node, target_comm_id, T=5):
        """
        Apply progressive BIH algorithm that hides from one community at a time
        while strengthening membership in the target community
        """
        G_modified = G.copy()
        target_community = communities[target_comm_id]
        other_communities = {cid: nodes for cid, nodes in communities.items() 
                             if cid != target_comm_id and target_node in nodes}
        
        operations = []
        
        # Calculate community metrics with more features
        comm_metrics = {}
        for comm_id, nodes in {target_comm_id: target_community, **other_communities}.items():
            subg = G.subgraph(nodes)
            size = len(nodes)
            density = nx.density(subg)
            edges = subg.number_of_edges()
            comm_metrics[comm_id] = {
                'size': size,
                'density': density,
                'avg_degree': 2 * edges / size if size > 0 else 0,
                'clustering': nx.average_clustering(subg) if size > 2 else 0
            }

        # Progressive hiding strategy:
        # 1. Order communities by overlap strength (most overlapping first)
        ordered_comms = []
        for other_comm, nodes in other_communities.items():
            overlap = len(set(G.neighbors(target_node)) & nodes)
            ordered_comms.append((overlap, other_comm))
        
        ordered_comms.sort(reverse=True)  # Handle strongest overlaps first
        
        # 2. Progressively hide from each overlapping community
        for _, other_comm in ordered_comms:
            # Hide from this community while strengthening target community
            G_temp, ops = self.hide_from_single_community(
                G_modified, 
                target_node,
                other_comm,
                target_comm_id,
                communities,
                comm_metrics
            )
            
            # If the change looks good (reduced overlap), keep it
            orig_overlap = len(set(G_modified.neighbors(target_node)) & communities[other_comm])
            new_overlap = len(set(G_temp.neighbors(target_node)) & communities[other_comm])
            
            if new_overlap < orig_overlap:
                G_modified = G_temp
                operations.extend(ops)
        
        # Enhanced edge removal phase with adaptive thresholds
        for comm_id, comm_nodes in other_communities.items():
            metrics = comm_metrics[comm_id]
            neighbors_in_comm = set(G_modified.neighbors(target_node)) & comm_nodes
            
            # Adaptive minimum edges based on community structure
            min_edges = max(1, int(metrics['density'] * len(neighbors_in_comm) * 0.3))
            
            while len(neighbors_in_comm) > min_edges:
                best_remove = max(
                    neighbors_in_comm,
                    key=lambda n: (
                        self.calculate_importance_degree(G_modified, n, comm_nodes) *
                        (1 + metrics['clustering']) * # Higher clustering = easier to remove
                        (G_modified.degree(n) / metrics['avg_degree'])
                    )
                )
                G_modified.remove_edge(target_node, best_remove)
                operations.append(f"Remove: {target_node}-{best_remove}")
                neighbors_in_comm.remove(best_remove)
        
        # Final phase: Add any remaining strategic edges to target community if needed
        target_metrics = comm_metrics[target_comm_id]
        current_neighbors = set(G_modified.neighbors(target_node))
        target_neighbors = current_neighbors & target_community
        
        if len(target_neighbors) < target_metrics['avg_degree'] * 0.7:
            potential_adds = target_community - current_neighbors - {target_node}
            potential_adds = sorted(
                potential_adds,
                key=lambda n: (
                    len(set(G_modified.neighbors(n)) & target_community) / (G_modified.degree(n) + 1),
                    -len(set(G_modified.neighbors(n)) & set().union(*[nodes for cid, nodes in other_communities.items()]))
                ),
                reverse=True
            )
            
            needed = int(target_metrics['avg_degree'] * 0.7) - len(target_neighbors)
            for n in potential_adds[:needed]:
                G_modified.add_edge(target_node, n)
                operations.append(f"Add: {target_node}-{n}")
        
        return G_modified, operations

    def select_best_target_and_apply(self, G, communities, target_node):
        """
        Try each candidate target community for the overlapping node by
        simulating the BIH operations and selecting the one that yields the
        smallest estimated number of remaining memberships for the node.

        This is a lightweight simulation (uses neighborhood-based membership
        estimate) and keeps the best simulated graph and operations.
        """
        candidate_comms = [cid for cid, nodes in communities.items() if target_node in nodes]
        if not candidate_comms:
            return G, [], None, []

        best_score = float('inf')
        best_G = G
        best_ops = []
        best_comm = None
        candidate_scores = []

        # Order candidates by overlap strength with target node (most overlapping first)
        cand_scores = []
        for cid in candidate_comms:
            overlap = len(set(G.neighbors(target_node)) & communities[cid])
            cand_scores.append((overlap, cid))
        cand_scores.sort(reverse=True)
        ordered_cids = [cid for _, cid in cand_scores]

        # Optionally limit number of candidate simulations to top-K
        if self.max_candidate_sims is not None:
            ordered_cids = ordered_cids[:self.max_candidate_sims]

        # We'll try single-community simulations and also simple two-step pairs
        pair_scores = []
        for cid in ordered_cids:
            # Simulate BIH for this target community
            G_sim, ops = self.apply_bih_to_snapshot(G, communities, target_node, cid, T=3)
            # Better estimation: run a lightweight community detection on the
            # local neighborhood after rewiring and count how many local
            # communities include the target node.
            try:
                local_nodes = set(G_sim.neighbors(target_node)) | {target_node}
                if len(local_nodes) > 0:
                    G_local = G_sim.subgraph(local_nodes).copy()
                    # Import here to avoid importing heavy modules at top-level in some environments
                    from networkx.algorithms.community import greedy_modularity_communities
                    comms_local = list(greedy_modularity_communities(G_local))
                    est_memberships = sum(1 for c in comms_local if target_node in c)
                else:
                    est_memberships = 0
            except Exception:
                # Fallback to simple neighborhood overlap count
                new_neighbors = set(G_sim.neighbors(target_node))
                est_memberships = sum(1 for cc, nodes in communities.items() if len(new_neighbors & nodes) > 0)

            # Record candidate score
            candidate_scores.append({
                'candidate_comm': cid,
                'estimated_memberships': int(est_memberships),
                'operations': len(ops),
                'type': 'single'
            })

            # Prefer fewer memberships, then fewer operations
            score = est_memberships
            if score < best_score or (score == best_score and len(ops) < len(best_ops)):
                best_score = score
                best_G = G_sim
                best_ops = ops
                best_comm = cid

        # Now try simple 2-step pairs (apply BIH to cid1 then to cid2)
        # Only try pairs from the limited ordered_cids to keep runtime reasonable
        import itertools
        for cid1, cid2 in itertools.permutations(ordered_cids, 2):
            # simulate first step
            G_step1, ops1 = self.apply_bih_to_snapshot(G, communities, target_node, cid1, T=2)
            # simulate second step on result (note: communities argument still refers to original partitions)
            G_step2, ops2 = self.apply_bih_to_snapshot(G_step1, communities, target_node, cid2, T=2)

            try:
                local_nodes = set(G_step2.neighbors(target_node)) | {target_node}
                if len(local_nodes) > 0:
                    G_local = G_step2.subgraph(local_nodes).copy()
                    from networkx.algorithms.community import greedy_modularity_communities
                    comms_local = list(greedy_modularity_communities(G_local))
                    est_memberships = sum(1 for c in comms_local if target_node in c)
                else:
                    est_memberships = 0
            except Exception:
                new_neighbors = set(G_step2.neighbors(target_node))
                est_memberships = sum(1 for cc, nodes in communities.items() if len(new_neighbors & nodes) > 0)

            total_ops = len(ops1) + len(ops2)
            pair_scores.append({
                'candidate_pair': (cid1, cid2),
                'estimated_memberships': int(est_memberships),
                'operations': total_ops
            })

            # Consider pair if better than current best
            if est_memberships < best_score or (est_memberships == best_score and total_ops < len(best_ops)):
                best_score = est_memberships
                best_G = G_step2
                best_ops = ops1 + ops2
                best_comm = (cid1, cid2)

        # attach pair scores to candidate_scores for debug
        for p in pair_scores:
            candidate_scores.append({
                'candidate_comm': p.get('candidate_pair'),
                'estimated_memberships': p.get('estimated_memberships'),
                'operations': p.get('operations'),
                'type': 'pair'
            })

        return best_G, best_ops, best_comm, candidate_scores
    
    def save_modified_network(self, G, snapshot_id, output_file):
        """Save modified network in ncol format"""
        with open(output_file, 'a') as f:
            for edge in G.edges():
                f.write(f"{edge[0]}\t{edge[1]}\t{snapshot_id}\n")
    
    def run_archangel(self, input_file, output_path, prefix="original_detection"):
        """Run ArchAngel on a network file"""
        print(f"Running ArchAngel...")
        
        # Record start time so we can detect files created by this ArchAngel run
        start_ts = time.time()

        aa = a.ArchAngel(
            network_filename=input_file,
            threshold=self.threshold,
            match_threshold=self.match_threshold,
            min_comsize=self.min_comsize,
            save=True,
            outfile_path=output_path
        )
        aa.execute()
        print(f"Analysis complete!")
        
        # Find generated files containing 'ArchAngel_coms_' either in the output path
        # or in the current working directory (ArchAngel may write files without
        # a separating path). Collect full paths so we can rename them reliably.
        candidates = []  # list of (folder, filename)
        known_prefixes = ('original_detection', 'modified_detection')
        for folder in (output_path, os.getcwd()):
            try:
                for f in os.listdir(folder):
                    full = os.path.join(folder, f)
                    # Only consider text files that include the ArchAngel marker
                    if not (f.endswith('.txt') and 'ArchAngel_coms_' in f):
                        continue
                    # Skip files that are already prefixed by known prefixes
                    if any(f.startswith(p) for p in known_prefixes):
                        continue
                    # Use modification time to ensure file was created/updated by this run
                    try:
                        mtime = os.path.getmtime(full)
                    except OSError:
                        continue
                    if mtime >= start_ts - 1.0:  # allow 1s clock skew
                        candidates.append((folder, f))
            except FileNotFoundError:
                continue

        renamed_files = []
        for folder, filename in candidates:
            old_path = os.path.join(folder, filename)
            # Build new filename by replacing the leading chunk up to 'ArchAngel_coms_' with the prefix
            idx = filename.find('ArchAngel_coms_')
            if idx == -1:
                continue
            new_filename = prefix + filename[idx:]
            new_path = os.path.join(output_path, new_filename)

            # Ensure destination directory exists
            os.makedirs(output_path, exist_ok=True)

            # Remove existing target if present
            if os.path.exists(new_path):
                os.remove(new_path)

            # Move/rename file into the output_path with the new name
            os.rename(old_path, new_path)
            renamed_files.append(new_filename)

        # Also rename the matches CSV file if present (same replacement logic)
        for folder in (output_path, os.getcwd()):
            try:
                for f in os.listdir(folder):
                    if f.endswith('ArchAngel_coms_ct_matches.csv'):
                        old_matches = os.path.join(folder, f)
                        idx = f.find('ArchAngel_coms_ct_matches.csv')
                        new_matches_name = prefix + f[idx:]
                        new_matches = os.path.join(output_path, new_matches_name)
                        if os.path.exists(new_matches):
                            os.remove(new_matches)
                        os.rename(old_matches, new_matches)
            except FileNotFoundError:
                continue

        return sorted(renamed_files), prefix
    
    def extract_snapshot_id_from_filename(self, filename, prefix):
        """Extract snapshot ID from filename"""
        try:
            # Remove prefix and suffix
            snapshot_str = filename.replace(f'{prefix}ArchAngel_coms_', '').replace('.txt', '')
            return int(snapshot_str)
        except:
            return None
    
    def compare_results(self, original_comm_file, modified_comm_file, 
                          overlapping_nodes, snapshot_id, G_original):
        """Compare community structures before and after BIH"""
        print(f"\n   Loading communities for comparison:")
        print(f"       Original: {original_comm_file}")
        print(f"       Modified: {modified_comm_file}")
        
        orig_comms = self.load_communities(original_comm_file)
        mod_comms = self.load_communities(modified_comm_file)
        
        print(f"       Original communities loaded: {len(orig_comms)}")
        print(f"       Modified communities loaded: {len(mod_comms)}")
        
        results = []
        
        for node, orig_comm_ids in overlapping_nodes.items():
            print(f"\n         Comparing node {node}:")
            # Count communities before
            orig_count = len(orig_comm_ids)
            print(f"           Original communities ({orig_count}): {orig_comm_ids}")
            
            # Count communities after
            mod_comm_ids = [cid for cid, nodes in mod_comms.items() if node in nodes]
            mod_count = len(mod_comm_ids)
            print(f"           Modified communities ({mod_count}): {mod_comm_ids}")
            
            # Calculate hiding success
            hiding_success = mod_count == 1
            print(f"           Hiding {'successful' if hiding_success else 'failed'}")
            
            # Get degree
            degree = G_original.degree(node) if node in G_original else 0
            
            results.append({
                'snapshot': snapshot_id,
                'node': node,
                'degree_before': degree,
                'communities_before': orig_count,
                'communities_after': mod_count,
                'hiding_successful': hiding_success,
                'original_communities': ','.join(map(str, sorted(orig_comm_ids))),
                'modified_communities': ','.join(map(str, sorted(mod_comm_ids)))
            })
        
        return results
    
    def process_all_snapshots(self, input_file):
        """Main processing pipeline"""
        print(f"\n{'='*70}")
        print("BIH DYNAMIC NETWORK HIDING PIPELINE")
        print(f"{'='*70}\n")
        
        # Determine number of snapshots
        snapshots = set()
        with open(input_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    snapshots.add(int(parts[2]))
        
        print(f"Found {len(snapshots)} snapshots: {sorted(snapshots)}\n")
        
        # Step 1: Run ArchAngel on original network
        print("="*70)
        print("STEP 1: Analyzing Original Network")
        print("="*70)
        orig_comm_files, orig_prefix = self.run_archangel(
            input_file, self.output_dir, prefix="original_detection"
        )
        print(f"Generated {len(orig_comm_files)} community files\n")
        
        # Step 2: Process each snapshot with BIH
        print("="*70)
        print("STEP 2: Applying BIH to Each Snapshot")
        print("="*70)
        
        bih_modified_file = os.path.join(self.networks_dir, "merged_snapshots_modified.ncol")
        if os.path.exists(bih_modified_file):
            os.remove(bih_modified_file)
        
        processed_snapshots = {}
        total_operations = 0
        
        # Process each community file from original analysis
        for comm_filename in orig_comm_files:
            snapshot_id = self.extract_snapshot_id_from_filename(comm_filename, orig_prefix)
            if snapshot_id is None:
                continue
            
            print(f"\n{'─'*70}")
            print(f"Snapshot {snapshot_id}")
            print(f"{'─'*70}")
            
            orig_comm_file = os.path.join(self.output_dir, comm_filename)
            communities = self.load_communities(orig_comm_file)
            overlapping = self.find_overlapping_nodes(communities)
            
            print(f"   Communities: {len(communities)}")
            print(f"   Overlapping nodes: {len(overlapping)}")
            
            # Load the network for this snapshot
            G = self.load_snapshot_network(input_file, snapshot_id)
            print(f"   Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
            
            if not overlapping:
                print("   No overlapping nodes - copying original network")
                self.save_modified_network(G, snapshot_id, bih_modified_file)
                processed_snapshots[snapshot_id] = {
                    'overlapping': overlapping,
                    'target_node': None,
                    'G_original': G
                }
                continue
            
            # Find node with highest degree among overlapping nodes
            degrees = {node: G.degree(node) for node in overlapping.keys() if node in G}
            
            if not degrees:
                print("   Warning: No overlapping nodes found in graph")
                self.save_modified_network(G, snapshot_id, bih_modified_file)
                processed_snapshots[snapshot_id] = {
                    'overlapping': overlapping,
                    'target_node': None,
                    'G_original': G
                }
                continue
            
            target_node = max(degrees, key=degrees.get)
            target_degree = degrees[target_node]
            target_comms = overlapping[target_node]

            print(f"   Target node: {target_node}")
            print(f"     Degree: {target_degree}")
            print(f"     Current communities: {target_comms}")

            # Improved selection: simulate BIH for each candidate target community
            if self.simulate_targets:
                G_modified, operations, chosen_comm, candidate_scores = self.select_best_target_and_apply(
                    G, communities, target_node
                )
            else:
                # Fallback: pick first community and run BIH once
                chosen_comm = target_comms[0]
                G_modified, operations = self.apply_bih_to_snapshot(
                    G, communities, target_node, chosen_comm, T=3
                )

            if chosen_comm is None:
                # Fallback to first community
                chosen_comm = target_comms[0]

            print(f"     Chosen target community: {chosen_comm}")

            # If debug enabled and we have candidate scores, save them for inspection
            if getattr(self, 'debug_sim', False) and self.simulate_targets:
                try:
                    debug_path = os.path.join(self.comparison_dir, f"sim_debug_snapshot_{snapshot_id}.json")
                    with open(debug_path, 'w') as jd:
                        json.dump({'snapshot': snapshot_id, 'target_node': target_node, 'candidates': candidate_scores}, jd, indent=2)
                except Exception:
                    pass
            
            if operations:
                print(f"   Applied {len(operations)} operations:")
                for op in operations[:5]:  # Show first 5 operations
                    print(f"         • {op}")
                if len(operations) > 5:
                    print(f"         ... and {len(operations)-5} more")
                total_operations += len(operations)
            else:
                print("   No operations needed")
            
            # Save modified network
            self.save_modified_network(G_modified, snapshot_id, bih_modified_file)
            
            processed_snapshots[snapshot_id] = {
                'overlapping': overlapping,
                'target_node': target_node,
                'G_original': G,
                'operations': len(operations)
            }
        
        print(f"\n{'─'*70}")
        print(f"Modified network saved: {bih_modified_file}")
        print(f"Total BIH operations across all snapshots: {total_operations}\n")
        
        # Step 3: Run ArchAngel on modified network
        print("="*70)
        print("STEP 3: Analyzing BIH-Modified Network")
        print("="*70)
        mod_comm_files, mod_prefix = self.run_archangel(
            bih_modified_file, self.output_dir, prefix="modified_detection"
        )
        print(f"Generated {len(mod_comm_files)} community files\n")
        
        # Step 4: Compare results
        print("="*70)
        print("STEP 4: Comparing Results")
        print("="*70)
        
        all_comparison_results = []
        
        for comm_filename in orig_comm_files:
            snapshot_id = self.extract_snapshot_id_from_filename(comm_filename, orig_prefix)
            if snapshot_id is None or snapshot_id not in processed_snapshots:
                continue
            
            orig_comm_file = os.path.join(self.output_dir, comm_filename)
            mod_comm_filename = comm_filename.replace(orig_prefix, mod_prefix)
            mod_comm_file = os.path.join(self.output_dir, mod_comm_filename)
            
            print(f"\nSnapshot {snapshot_id} comparison:")
            print(f"   Original file: {orig_comm_file}")
            print(f"   Modified file: {mod_comm_file}")

            if not os.path.exists(mod_comm_file):
                print(f"Warning: Modified community file not found for snapshot {snapshot_id}")
                continue
            
            snapshot_data = processed_snapshots[snapshot_id]
            overlapping = snapshot_data['overlapping']
            G_original = snapshot_data['G_original']
            target_node = snapshot_data.get('target_node')
            
            print(f"   Overlapping nodes: {len(overlapping)}")
            print(f"   Target node: {target_node}")
            if target_node:
                print(f"   Target communities: {overlapping.get(target_node, [])}")
            
            # Only compare the single targeted overlapping node (highest-degree)
            if overlapping and target_node is not None and target_node in overlapping:
                filtered_overlapping = {target_node: overlapping[target_node]}
                print(f"   Comparing target node {target_node}")
                results = self.compare_results(
                    orig_comm_file, mod_comm_file, filtered_overlapping, snapshot_id, G_original
                )
                print(f"   Comparison results: {len(results)} rows")
                for r in results:
                    print(f"       {r['node']}: {r['communities_before']} → {r['communities_after']} communities")
                    # Extend results with snapshot data
                    # *** Add the number of operations to the result dictionary ***
                    r['operations'] = processed_snapshots[snapshot_id].get('operations', 0)
                    all_comparison_results.extend(results)
            else:
                print(f"   Skipping comparison - conditions not met:")
                print(f"     Has overlapping: {bool(overlapping)}")
                print(f"     Has target_node: {target_node is not None}")
                print(f"     Target in overlapping: {target_node in overlapping if target_node else False}")
                results = []  # Print results for this snapshot
                print(f"\nSnapshot {snapshot_id}:")
                for result in results:
                    status = "SUCCESS" if result['hiding_successful'] else "FAILED"
                    print(f"   Node {result['node']}: "
                          f"{result['communities_before']} → {result['communities_after']} "
                          f"communities ({status})")
                
                all_comparison_results.extend(results)
        
        # Save comparison results
        comparison_file = os.path.join(self.comparison_dir, "hiding_comparison.csv")
        
        # Ensure we write results even if some snapshots had no comparisons
        if len(all_comparison_results) > 0:
            # Get fieldnames from the first result, ensuring 'operations' is included
            fieldnames = list(all_comparison_results[0].keys())
            if 'operations' not in fieldnames:
                fieldnames.append('operations')

            with open(comparison_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_comparison_results)
            
            # Print detailed summary
            print("\n" + "="*70)
            print("FINAL SUMMARY")
            print("="*70)
            
            successful = sum(1 for r in all_comparison_results if r['hiding_successful'])
            total = len(all_comparison_results)
            
            print(f"\nHiding Results:")
            print(f"   Total target nodes processed: {total}")
            print(f"   Successfully hidden (moved to single community): {successful}")
            print(f"   Failed to hide: {total - successful}")
            print(f"   Success rate: {successful/total*100:.1f}%")
            
            # Detailed breakdown
            by_snapshot = defaultdict(list)
            for result in all_comparison_results:
                by_snapshot[result['snapshot']].append(result)
            
            print(f"\nBreakdown by Snapshot:")
            for snapshot in sorted(by_snapshot.keys()):
                results = by_snapshot[snapshot]
                success = sum(1 for r in results if r['hiding_successful'])
                print(f"   Snapshot {snapshot}: {success}/{len(results)} successful")
        else:
            # Write empty CSV with headers
            with open(comparison_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['snapshot', 'node', 'degree_before', 
                                                       'communities_before', 'communities_after',
                                                       'hiding_successful', 'original_communities',
                                                       'modified_communities', 'operations'])
                writer.writeheader()
            
            print("\n" + "="*70)
            print("FINAL SUMMARY")
            print("="*70)
            print("\nNo target nodes were compared")
        
        print(f"\nOutput Files:")
        print(f"   Original communities: {self.output_dir}")
        print(f"   Modified network: {bih_modified_file}")
        print(f"   Comparison CSV: {comparison_file}")

        # --- NEW: Call the plotting function ---
        if PLOTLIBS_AVAILABLE and len(all_comparison_results) > 0:
            self.plot_hiding_cost_per_snapshot(all_comparison_results)
        elif PLOTLIBS_AVAILABLE:
            print("\nPlotting skipped: No results to plot.")
        else:
            print("\nPlotting skipped: 'pandas' or 'matplotlib' not installed.")
        # --- End of new section ---

        print("="*70 + "\n")


    # --- NEW FUNCTION ---
    def plot_hiding_cost_per_snapshot(self, all_comparison_results):
        """
        Generates and saves a bar chart of the hiding cost (operations) per snapshot.
        """
        if not PLOTLIBS_AVAILABLE:
            print("Error: Plotting libraries (pandas, matplotlib) not available.")
            return

        try:
            df = pd.DataFrame(all_comparison_results)
            
            if 'snapshot' not in df.columns or 'operations' not in df.columns:
                print("Warning: Could not plot hiding cost. 'snapshot' or 'operations' column missing.")
                return

            # Ensure data is sorted by snapshot for a clean plot
            df = df.sort_values(by='snapshot')

            # Get unique snapshots and corresponding operations
            plot_data = df.drop_duplicates(subset=['snapshot'])

            plt.figure(figsize=(10, 6))
            plt.bar(plot_data['snapshot'], plot_data['operations'], color='c', edgecolor='black')
            
            plt.xlabel('Snapshot ID')
            plt.ylabel('Number of Operations (Edge Adds/Removes)')
            plt.title('Hiding Cost (Operations) per Snapshot')
            plt.xticks(plot_data['snapshot'].astype(int))
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            
            # Save the figure
            save_path = os.path.join(self.output_dir, "hiding_cost_per_snapshot.png")
            plt.savefig(save_path)
            plt.close()
            
            print(f"Hiding cost chart saved to: {save_path}")

        except Exception as e:
            print(f"Error: Could not generate hiding cost plot. {e}")
    # --- END OF NEW FUNCTION ---


# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run BIH dynamic hiding pipeline')
    parser.add_argument('input_file',
                        help='Path to the input .ncol file containing the dynamic network')
    parser.add_argument('--base-path', default=os.getcwd(),
                        help='Base path for output folders (default: current directory)')
    parser.add_argument('--no-simulate-targets', dest='simulate_targets', action='store_false',
                        help='Disable simulation-based target selection')
    parser.add_argument('--debug-sim', dest='debug_sim', action='store_true',
                        help='Write per-snapshot simulation debug files')
    parser.add_argument('--max-candidate-sims', dest='max_candidate_sims', type=int, default=None,
                        help='Limit number of candidate communities to simulate (top-K)')
    parser.add_argument('--no-hillclimb', dest='enable_hillclimb', action='store_false',
                        help='Disable hill-climb local search')
    parser.add_argument('--hillclimb-max-iter', dest='hillclimb_max_iter', type=int, default=3,
                        help='Max iterations for hill-climb local search')
    args = parser.parse_args()

    base_path = args.base_path

    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}")
        sys.exit(1)

    # Initialize and run BIH pipeline
    bih = BIHDynamicHiding(base_path,
                           simulate_targets=args.simulate_targets,
                           debug_sim=args.debug_sim,
                           max_candidate_sims=args.max_candidate_sims,
                           enable_hillclimb=args.enable_hillclimb,
                           hillclimb_max_iter=args.hillclimb_max_iter)
    bih.process_all_snapshots(args.input_file)