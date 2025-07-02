import dendropy

import csv
import re
import json
import argparse
from typing import FrozenSet
import os

# Define a safe placeholder for commas within 'rea' annotations
COMMA_PLACEHOLDER = "||COMMA||"


def fix_and_parse_tree(tree_path, is_original=False):
    """
    Reads a tree file, performs a safe replacement of commas inside 'rea'
    annotations, and then parses the fixed string with dendropy.
    If is_original, it also labels unlabeled internal nodes.
    """
    if not os.path.exists(tree_path):
        print(f"Error: Input file not found at '{tree_path}'")
        return None

    with open(tree_path, "r") as f:
        tree_string = f.read()

    def replace_commas_in_match(match: re.Match) -> str:
        content_with_placeholders = match.group(1).replace(",", COMMA_PLACEHOLDER)
        return f'rea="{content_with_placeholders}"'

    # This is a regular expression to replace commas in Treesort node annotations with a placeholder
    fixed_tree_string = re.sub(r"rea=\"(.*?)\"", replace_commas_in_match, tree_string)

    schema = "newick" if is_original else "nexus"
    tree = dendropy.Tree.get(
        data=fixed_tree_string, schema=schema, preserve_underscores=True
    )

    if is_original:
        counter = 0
        # Iterate in post-order to give child nodes lower numbers than parents
        for node in tree.postorder_internal_node_iter():
            if node.label is None:
                node.label = f"NODE_{counter}"
                counter += 1
    return tree


def get_rea_string_from_node(node):
    """
    Finds the 'rea' annotation string for a node by checking all likely locations,
    and replaces the placeholder back to a comma.
    """
    locations_to_check = []
    if node.edge:
        locations_to_check.append(node.edge.annotations)
    locations_to_check.append(node.annotations)
    if node.taxon:
        locations_to_check.append(node.taxon.annotations)

    for annotation_set in locations_to_check:
        for item in annotation_set:
            if item.name.lower() == "rea":
                return item.value.replace(COMMA_PLACEHOLDER, ",")
    return None


def get_leaf_set_map(tree):
    """
    Builds a dictionary mapping a node's label to the frozenset of its leaves.
    This works for BOTH internal and leaf nodes.
    """
    clade_map = {}
    for node in tree.postorder_node_iter():
        node_label = node.label if node.label else node.taxon.label
        if node_label:
            leaf_labels = {leaf.taxon.label for leaf in node.leaf_nodes()}
            clade_map[node_label] = frozenset(leaf_labels)
    return clade_map


def write_merged_report_to_csv(report, csv_path):
    """
    Converts the hierarchical report into a final CSV file formatted for Augur.
    """
    print(f"\n--- Creating final merged CSV output at '{csv_path}' ---")

    header = ["strain", "reassortment_events", "source"]
    rows = []

    for node_label, data in report.items():
        raw_event_strings = []
        source_strings = []

        if data.get("Direct Annotations"):
            source_strings.append("DIRECT")
            raw_event_strings.extend(data["Direct Annotations"])

        if data.get("Subclade Annotations"):
            source_strings.append("SUBCLADE")
            for sub_event in data["Subclade Annotations"]:
                raw_event_strings.append(sub_event["event"])

        if raw_event_strings:
            final_events = ",".join(raw_event_strings)
            final_sources = " | ".join(
                sorted(list(set(source_strings)))
            )  # ensure consistent order

            rows.append(
                {
                    "strain": node_label,
                    "reassortment_events": final_events,
                    "source": final_sources,
                }
            )

    if rows:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Successfully wrote {len(rows)} unique node annotations to {csv_path}.")
    else:
        print("Warning: No mapped annotations to write to CSV.")


def run_pipeline(args):
    """Main function to run the simplified analysis pipeline."""

    # --- Part 1: Safely parse the trees and label the original ---
    processed_tree = fix_and_parse_tree(args.treesort_nexus)
    if not processed_tree:
        return

    original_tree = fix_and_parse_tree(args.original_tree, is_original=True)
    if not original_tree:
        return

    original_tree.write(path=args.labeled_tree, schema="newick", suppress_rooting=True)
    print(f"\nSuccessfully created labeled original tree at: {args.labeled_tree}")

    # --- Part 2: Perform Hierarchical Node Mapping ---
    print(f"\n--- Performing Hierarchical Node Mapping ---")

    # Create the leaf set maps for both trees
    processed_clades = get_leaf_set_map(processed_tree)
    original_clades = get_leaf_set_map(original_tree)

    # Invert original map for finding matches by set
    original_sets_to_labels = {v: k for k, v in original_clades.items()}

    report = {}

    # Iterate through all annotations found in the processed tree
    for proc_node_label, proc_leaf_set in processed_clades.items():
        # Find the annotation for this node
        rea_string = None
        proc_node = processed_tree.find_node_with_label(proc_node_label)
        if not proc_node:  # It might be a leaf node
            proc_node = processed_tree.find_node_with_taxon_label(proc_node_label)
        if proc_node:
            rea_string = get_rea_string_from_node(proc_node)

        if not rea_string:
            continue

        # --- Logic for ALL nodes (leaves and internal) ---

        # Scenario A: Direct Match
        if proc_leaf_set in original_sets_to_labels:
            original_node_label = original_sets_to_labels[proc_leaf_set]
            if original_node_label not in report:
                report[original_node_label] = {
                    "Direct Annotations": [],
                    "Subclade Annotations": [],
                }
            report[original_node_label]["Direct Annotations"].append(rea_string)

        # Scenario B: Subclade Match (for internal nodes only)
        elif len(proc_leaf_set) > 1:
            best_container_label = None
            min_size = float("inf")

            for orig_label, orig_set in original_clades.items():
                if proc_leaf_set.issubset(orig_set):
                    if len(orig_set) < min_size:
                        min_size = len(orig_set)
                        best_container_label = orig_label

            if best_container_label:
                if best_container_label not in report:
                    report[best_container_label] = {
                        "Direct Annotations": [],
                        "Subclade Annotations": [],
                    }
                report[best_container_label]["Subclade Annotations"].append(
                    {
                        "event": rea_string,
                        "on_subclade_defined_by_leaves": sorted(list(proc_leaf_set)),
                    }
                )

    print("Hierarchical mapping complete.")

    if report:
        write_merged_report_to_csv(report, args.treesort_csv)
    else:
        print("\nWarning: Final report was empty. No mapped CSV was generated.")

    # Generate JSON output
    if report:
        convert_csv_to_auspice_json(
            args.treesort_csv, args.labeled_tree, args.treesort_json
        )


def convert_csv_to_auspice_json(csv_path, tree_path, json_output_path):
    """
    Converts the mapped reassortment CSV into a JSON format for Auspice,
    containing both 'nodes' and 'branches' attributes.
    """
    print(f"--- Converting '{csv_path}' to Auspice JSON format ---")
    if not os.path.exists(csv_path):
        print(f"Error: CSV input file not found at '{csv_path}'")
        return
    if not os.path.exists(tree_path):
        print(f"Error: Tree file not found at '{tree_path}'")
        return

    reassortment_data = {}
    try:
        with open(csv_path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                reassortment_data[row["strain"]] = {
                    "reassortment_events": row["reassortment_events"],
                    "source": row["source"],
                }
    except Exception as e:
        print(f"An error occurred reading the CSV file: {e}")
        return

    try:
        # Parse the tree and extract all node names directly
        tree = dendropy.Tree.get_from_path(
            tree_path, schema="newick", preserve_underscores=True
        )
        all_node_names = []
        for node in tree.postorder_node_iter():
            if node.taxon and node.taxon.label:
                all_node_names.append(node.taxon.label)
            elif node.label:
                all_node_names.append(node.label)
    except Exception as e:
        print(f"An error occurred reading the tree file: {e}")
        return

    final_json = {
        "branch_attrs": {"labels": {"text": "Reassorted Segments"}},
        "nodes": {},
        "branches": {},
    }

    for node_name in all_node_names:
        if node_name in reassortment_data:
            final_json["nodes"][node_name] = {
                "Reassorted": "True",
                "reassortment_events": reassortment_data[node_name][
                    "reassortment_events"
                ],
                "source": reassortment_data[node_name]["source"],
            }

    for node_name, data in reassortment_data.items():
        # Extract segment names from patterns like "PB2(344),PB1(319)" or "NS(49)"
        # Split by comma, then extract the part before the opening parenthesis
        event_string = data["reassortment_events"]
        segments = []
        for part in event_string.split(","):
            part = part.strip().strip('"')  # Remove quotes and whitespace
            if "(" in part:
                segment = part.split("(")[0].lstrip("?_")  # Remove leading ? or _
                segments.append(segment)
        unique_segments = sorted(list(set(segments)))

        count = len(unique_segments)
        segments_str = ", ".join(unique_segments)

        final_json["branches"][node_name] = {
            "labels": {"Reassorted Segments": f"{count} ({segments_str})"}
        }

    with open(json_output_path, "w") as f:
        json.dump(final_json, f, indent=2)

    print(f"Successfully created Auspice JSON at: {json_output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A robust script to parse TreeSort results and map them to a newly-labeled original tree for Augur."
    )
    parser.add_argument("treesort_nexus", help="TreeSort-modified NEXUS file")
    parser.add_argument("original_tree", help="Original Newick tree file")
    parser.add_argument("labeled_tree", help="Output path for labeled original tree")
    parser.add_argument("treesort_csv", help="Output path for CSV file")
    parser.add_argument("treesort_json", help="Output path for Auspice JSON file")
    args = parser.parse_args()

    run_pipeline(args)
