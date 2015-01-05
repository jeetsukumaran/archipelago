#! /usr/bin/env python

try:
    import lagrange
    IS_LAGRANGE_AVAILABLE = True
except ImportError:
    IS_LAGRANGE_AVAILABLE = False

LAGRANGE_CONFIGURATION_TEMPLATE = """\
#!/usr/bin/env python
import os
import lagrange
data = \"\"\"\
### begin data
{{
 'area_adjacency': {area_adjacency},
 'area_dispersal': {area_dispersal},
 'area_labels': {area_labels},
 'base_rates': '__estimate__',
 'dispersal_durations': [10000.0],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20130526',
 'max_range_size': {max_range_size},
 'model_name': '{model_name}',
 'newick_trees': [{{'included': [],
                   'name': 'Tree0',
                   'newick': '{newick_tree_string}',
                   'root_age': {root_age}}}],
 'ranges': {ranges},
 'taxa': {taxon_name_list},
 'taxon_range_data': {taxon_range_data},
}}
### end data
\"\"\"

i = 0
while 1:
    if not i:
        outfname = "{model_name}.results.txt"
    else:
        outfname = "{model_name}.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)

"""

class LagrangeEstimator(object):

    def estimate_lagrange_rates(self,
            tree,
            profile_results):
        lagrange_commands = self.compose_lagrange_template(tree=tree)
        commandsf = open(self.commands_file_name, "w")
        commandsf.write(lagrange_commands)
        commandsf.flush()
        commandsf.close()
        shell_cmd = ["python", self.commands_file_name]
        p = subprocess.Popen(
                shell_cmd,
                stdout=subprocess.PIPE,
                )
        stdout, stderr = processio.communicate(p)

    def compose_lagrange_template(self, tree):
        area_names = sorted(self.reconstruct_areas(tree))
        kwargs = {}
        kwargs["area_adjacency"] = str([[1] * len(area_names)] * len(area_names))
        kwargs["area_dispersal"] = str([[1.0] * len(area_names)] * len(area_names))
        kwargs["area_labels"] = str(area_names)
        kwargs["max_range_size"] = len(area_names)
        kwargs["model_name"] = str(id(tree))
        kwargs["newick_tree_string"] = tree.as_string("newick").replace("\n", "")
        assert tree.seed_node.age
        kwargs["root_age"] = tree.seed_node.age

        kwargs["taxon_name_list"] = [taxon.label for taxon in tree.taxon_namespace]

        taxon_range_data = {}
        for taxon in tree.taxon_namespace:
            taxon_area_indexes = tuple([area_idx for area_idx, i in enumerate(taxon.distribution_vector) if str(i) == "1"])
            taxon_range_data[taxon.label] = taxon_area_indexes
        kwargs["taxon_range_data"] = taxon_range_data

        ## EVERY PERMUTATION OF AREAS
        # ranges = []
        # for i in range(len(area_names)):
        #     x = list(itertools.permutations(area_indexes, i))
        #     ranges.extend(x)
        # ranges = sorted(set(ranges))
        # kwargs["ranges"] = str(ranges)

        ranges = set()
        for a in taxon_range_data.values():
            ranges.add(a)
        ranges = sorted(ranges)
        kwargs["ranges"] = str(ranges)

        return LAGRANGE_CONFIGURATION_TEMPLATE.format(**kwargs)

    def reconstruct_areas(self, tree):
        # cannot rely on generating model for number of
        # areas because we do not know if supplemental
        # areas are included in node data
        num_areas = len(tree.taxon_namespace[0].distribution_vector)
        area_names = ["a{}".format(i+1) for i in range(num_areas)]
        return area_names

