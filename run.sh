python mapper.py data/summary_reassortment.json data/output.nwk input.nwk labeled.nwk node.json --debug > debug-log.txt
augur export v2 --tree labeled.nwk --metadata h3nx_ha.csv --node-data node.json --auspice-config ./config/auspice_treesort.json --colors config/colors_h3nx.tsv --lat-longs config/lat_longs_h3nx.tsv --output h3nx_treesort_auspice.json
