#! /usr/bin/env python

import os
import collections
import argparse
import subprocess
import csv
import re
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.interop import rstats
from dendropy.utility import cli
import archipelago
from archipelago import utility

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map

def main():
    parser = argparse.ArgumentParser(
            formatter_class=cli.CustomFormatter,
            )
    parser.add_argument(
            "target_summary_stats",
            help="Path to summary statistic data to be classified.")
    parser.add_argument(
            "training_summary_stats",
            nargs="+",
            help="Path(s) to training summary statistic data.")
    npca_options = parser.add_argument_group("Number of Principle Component Axes Retained",
                (
                "Number of axes retained in the principal component step."
                " Exactly one of '--set-npca' or '--calculate-npca'"
                " must be specified."
                ))
    # npca_options = npca_options_parent_group.add_mutually_exclusive_group(required=True)
    npca_options.add_argument("--set-npca",
            default=None,
            metavar="n",
            help="Number of axes retained in the principal component analysis step.",
            )
    npca_options.add_argument("--calculate-npca",
            default=None,
            metavar="CRITERIA",
            choices=["all", "maximized-fit", "penalized-fit"],
            help=("\n".join([
                "<pre>Calculate number of axes to be retained in the" ,
                "principal component step based on one of the following",
                "criteria:",
                " - 'all'",
                "     Use maximum number of principal component ",
                "     axes available (i.e., number of summary ",
                "     statistics). ",
                " - 'maximized-fit'",
                "     Use number of principal component axes ",
                "     that maximizes the proportion of correct ",
                "     classifications when reapplied to the training",
                "     data. ",
                " - 'penalized-fit'",
                "     Use number of principal component axes that ",
                "     maximizes a score given by the (sum of the ",
                "     proportion of correct classifications when ",
                "     reapplied to the training data) minus (the ",
                "     number of axes weighted by a penalty factor). ",
                ])))
    npca_options.add_argument("--fit-penalty-factor",
            default=0.1,
            type=float,
            help="Penalty factor for penalized fitting (default: %(default)s).")
    nda_options = parser.add_argument_group("Number of Discriminant Analysis Axes Retained",
                (
                "Number of axes retained in the discriminant analysis step."
                ))
    nda_options.add_argument("--nda",
            default=None,
            metavar="n",
            help=(
                "Number of axes retained in the discriminant analysis step."
                " If not specified, defaults to one less than the number of"
                " candidate models (groups) in the data."
            ))
    extended_analysis_options = parser.add_argument_group("Extended Analysis Options")
    extended_analysis_options.add_argument(
            "--true-model",
            default=None,
            help=(
                "The label or name of the true model for the target data. This is typically"
                " only known if the target data was simulated as well, and you are interested"
                " in assessing the performance of the inference method or simulation regime."
                " If specified, will result in performance assessment statistics being generated."
                " These will be written to standard output after the primary results, or to the"
                " filepath specified by '--performance-output'."
                ))
    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    output_options.add_argument(
            "-o", "--primary-output-filepath",
            default=None,
            metavar="FILEPATH",
            help="Path to primary results output file.")
    output_options.add_argument(
            "--no-primary-output",
            action="store_true",
            default=False,
            help="Do not write primary output.")
    output_options.add_argument(
            "--performance-output-filepath",
            default=None,
            metavar="FILEPATH",
            help="Path to performance results output file (will only be written if '--true-model' is specified.")
    output_options.add_argument(
            "--no-performance-output",
            action="store_true",
            default=False,
            help="Do not write performance output.")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append instead of overwriting output file(s).")
    args = parser.parse_args()

    # extra fields
    extra_fields = parse_fieldname_and_value(args.labels)

    # true model performance
    if args.true_model is not None:
        true_model_name = args.true_model
        is_performance_assessed = True
    else:
        true_model_name = None
        is_performance_assessed = False

    training_summary_stats_paths = "c({})".format(",".join("'{}'".format(f) for f in args.training_summary_stats))

    if args.set_npca is None and args.calculate_npca is None:
        sys.exit("Must specify one of '--set-npca' or '--calculate-npca'")
    elif args.set_npca is not None and args.calculate_npca is not None:
        sys.exit("Must specify only one of '--set-npca' or '--calculate-npca'")
    elif args.set_npca is not None:
        try:
            npca = int(args.set_npca)
        except ValueError:
            sys.exit("Invalid integer specified for '--npca': '{}'".format(args.set_npca))
    else:
        npca = "'{}'".format(args.calculate_npca)

    if args.nda is None:
        nda = "NULL"
    else:
        try:
            npca = int(args.set_npca)
        except ValueError:
            sys.exit("Invalid integer specified for '--npca': '{}'".format(args.set_npca))


    r_commands = []
    r_commands.append("source('{}')".format(archipelago.libexec_filepath("analyze-dapc.R")))
    r_commands.append("""
            results = classifyDataFromFiles(
                    target.summary.stats.path='{target_summary_stats_path}',
                    training.summary.stats.paths={training_summary_stats_paths},
                    n.pca={npca},
                    n.da={nda},
                    n.pca.optimization.penalty.weight={fit_penalty_factor},
                    output.path="")
            """.format(
                target_summary_stats_path=args.target_summary_stats,
                training_summary_stats_paths=training_summary_stats_paths,
                npca=npca,
                nda=nda,
                fit_penalty_factor=args.fit_penalty_factor,
                ))

    returncode, stdout, stderr = rstats.call(r_commands)
    print(stdout)
    # assert os.path.exists(r_script_path)
    # cmd = [
    #         r_script_path,
    #         args.npca,
    #         args.nda,
    #         args.target_summary_stats,]
    # cmd.extend(args.training_summary_stats)
    # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = p.communicate()
    # if p.returncode != 0:
    #     sys.exit(stderr)

    reader = csv.DictReader(StringIO(stdout))

    correct_assignment_count = 0.0
    sum_of_true_model_pp = 0.0
    total_assignments = 0
    primary_result_rows = []
    # correct_assignment_count_by_model = collections.defaultdict(lambda: 0.0)
    # incorrect_assignment_count_by_model = collections.defaultdict(lambda: 0.0)
    for row in reader:
        row.update(extra_fields)
        primary_result_rows.append(row)
        if is_performance_assessed:
            total_assignments += 1
            row["true.model"] = true_model_name
            row["true.model.posterior"] = row["posterior.{}".format(true_model_name)]
            sum_of_true_model_pp += float(row["true.model.posterior"])
            assigned_model_name = row["assign"]
            if assigned_model_name == true_model_name:
                row["is.correctly.assigned"] = "T"
                correct_assignment_count += 1
                # correct_assignment_count_by_model[true_model_name] += 1
            else:
                row["is.correctly.assigned"] = "F"
                # incorrect_assignment_count_by_model[true_model_name] += 1

    if not args.no_primary_output:
        primary_result_fieldnames = []
        primary_result_fieldnames.extend([f for f in extra_fields.keys() if f not in primary_result_fieldnames])
        primary_result_fieldnames.extend([f for f in reader.fieldnames if f not in primary_result_fieldnames])
        if is_performance_assessed:
            primary_result_fieldnames.extend([
                "true.model",
                "true.model.posterior",
                "is.correctly.assigned",
                ])
        out = utility.open_output_file_for_csv_writer(
                filepath=args.primary_output_filepath,
                append=args.append)
        writer = csv.DictWriter(out,
                fieldnames=primary_result_fieldnames,
                lineterminator=os.linesep,
                )
        if not args.append:
            writer.writeheader()
        writer.writerows(primary_result_rows)

    if is_performance_assessed and not args.no_performance_output:
        performance_row = collections.OrderedDict()
        performance_row.update(extra_fields)
        # performance_row["true.model"] = true_model_name
        performance_row["true.model.posterior.mean"] = sum_of_true_model_pp / total_assignments
        performance_row["true.model.proportion.correctly.assigned"] = correct_assignment_count / total_assignments
        # model_names = sorted(set( set(correct_assignment_count_by_model.keys()) | set(incorrect_assignment_count_by_model.keys()) ))
        # for model_name in model_names:
        #     performance_row["{}.proportion.correctly.assigned".format(model_name)] = correct_assignment_count_by_model[model_name]/total_assignments
        #     performance_row["{}.proportion.incorrectly.assigned".format(model_name)] = incorrect_assignment_count_by_model[model_name]/total_assignments
        if (
                (args.primary_output_filepath is None or args.primary_output_filepath == "-")
                and (args.performance_output_filepath is None or args.performance_output_filepath == "-")
                ):
            sys.stdout.write(chr(28)) # chr(28) = FS = File Separator
        out = utility.open_output_file_for_csv_writer(
                filepath=args.performance_output_filepath,
                append=args.append)
        writer = csv.DictWriter(out,
                fieldnames=performance_row.keys(),
                lineterminator=os.linesep,
                )
        if not args.append:
            writer.writeheader()
        writer.writerow(performance_row)

if __name__ == "__main__":
    main()
