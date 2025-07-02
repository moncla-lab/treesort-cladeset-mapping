# TreeSort Cladeset Mapping

This repository contains tools for mapping reassortment annotations from phylogenetic trees that have been modified by TreeSort back to the original unmodified trees. The purpose is support Nextstrain Auspice visualizations using TreeSort to identify reassortment events in influenza virus sequences.

## Overview

TreeSort is a tool that analyzes phylogenetic trees and multiple sequence alignments to identify reassortment events in segmented viruses like influenza. However, it can modify the tree by either binarizing polytomies or collapsing internal nodes. This can render trees difficult to view, particularly with Nextstrain and Auspice. There is thus a need for annotations to be mapped back to the original tree structure.

## Installation

While this tool will find most utility when incorporated with other repos, for reproducbility we provide enough context to run a standalone example.

Make sure you have installed [Bioconda](https://bioconda.github.io/) and configured it correctly as described in the link.

To install the environment:

```
conda create -n treesort-cladesets dendropy augur
```

Before running scripts, open a new shell and activate the environment:

```
conda activate treesort-cladesets
```

## Usage

**`mapper.py`** - Main processing script that:

- Parses TreeSort NEXUS files with reassortment annotations
- Maps events from modified trees back to original trees using cladeset matching
- Handles both direct matches (identical leaf sets) and subclade matches
- Outputs CSV files compatible with Augur for Nextstrain visualization
- Generates Auspice JSON format for direct visualization

```bash
python mapper.py treesort.nex input.nwk labeled_tree.nwk output.csv output.json
```

**Arguments (in order):**

1. `treesort_nexus` - TreeSort output NEXUS file
2. `original_tree` - Original Newick tree file  
3. `labeled_tree` - Output path for labeled original tree
4. `treesort_csv` - Output path for CSV file
5. `treesort_json` - Output path for Auspice JSON file

An example of how to visualize this with Nextstrain:

```
augur export v2 --tree input-labeled.nwk \
	--metadata h3nx_ha.csv \
	--node-data treesort-node-data.json \
	--auspice-config ./config/auspice_treesort.json \
	--colors config/colors_h3nx.tsv \
	--lat-longs config/lat_longs_h3nx.tsv \
	--output h3nx_treesort_auspice.json
```

## Details

The original tree is parsed, and every internal node is given a name to ensure consistency in the mapping step and the Auspice visualization. This tree should also be fed into Augur export to ensure that mapped reassortment events. Modifications in topology can be accommodated by cladeset matching.

### Cladeset matching

The purpose of this algorithm is to address the following: given a node in the TreeSort tree, how do we find what node should be assigned in the original tree, specifically for the tricky cases where the tree structures don't perfectly align? The following strategy utilized in this repo is built from the principle that a node is defined by the set of all tips that descend from it.

#### Step 1: Identify the Target
The process begins with a single reassortment event found on a node in the TreeSort tree. The script first determines the complete set of all leaves (or tips) that descend from this event node in the Treesort tree. This is our Target Set.

#### Step 2: Prepare for the Search
To find the best fit in the original tree, the script needs to keep track of the best match it has found so far. It creates two temporary placeholders:

One to store the label of the best-matching node (initially empty).

One to store the size (the number of leaves) of that best-matching node's clade (initially set to an infinitely large number).

#### Step 3: Loop Through All Possibilities (The Brute-Force Approach)
Now, the brute-force search begins. The script iterates through every single node in the original tree, one by one, to see if it's a potential match.

#### Step 4: The Containment Check
For each node in the original tree, the script gets its leaf set (we'll call this the Candidate Set). It then performs a simple but critical check: Does this Candidate Set contain all of the leaves from our Target Set?

If No, this original node is not a valid container for the reassortment event. The script immediately discards it and moves on to the next node in its loop.

If Yes, the node is a valid candidate, and the script proceeds to the next step to see if it's the best candidate so far.

#### Step 5: Find the Smallest Containing Clade (The "Tightest Fit")
Once the script has a valid candidate, it checks if it's a better fit than any previous match. It does this by comparing the size of the current Candidate Set to the size of the best match it has stored so far.

If the current candidate's leaf set is smaller than the best-so-far, it means we have found a "tighter" fit—a more recent common ancestor in the original tree that still contains all the necessary descendants. The script then updates its placeholders, storing the label and smaller size of this new best-matching node.

If the candidate's leaf set is larger than the best-so-far, it is ignored, as we have already found a more specific location for the event.

#### Step 6: Declare the Final Mapping
This process of checking every node in the original tree continues until all possibilities have been examined. At the end of the loop, the placeholder will hold the label of the single node in the original tree that represents the smallest possible clade that fully contains the Target Set. This node is declared the definitive mapping location for that reassortment event.

This exhaustive process guarantees that each event is mapped to the most specific and appropriate ancestral node in the original tree.