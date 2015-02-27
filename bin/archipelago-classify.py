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
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "target_summary_stats",
            help="Path to summary statistic data to be classified.")
    parser.add_argument(
            "training_summary_stats",
            nargs="+",
            help="Path(s) to training summary statistic data.")
    parser.add_argument(
            "--npca",
            default="NULL",
            help="Number of principle component axes to use in construction of the DAPC functions.")
    parser.add_argument(
            "--nda",
            default="NULL",
            help="Number of discriminant analysis functions to use.")
    parser.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    parser.add_argument(
            "-o", "--output-filepath",
            default=None,
            help="Path to output file.")
    parser.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append to output file if it already exists instead of overwriting.")
    args = parser.parse_args()
    extra_fields = parse_fieldname_and_value(args.labels)
    r_script_path = os.path.join(archipelago.ARCHIPELAGO_LIBEXEC_PATH, "archipelago-classify.R")
    assert os.path.exists(r_script_path)
    cmd = [
            r_script_path,
            args.npca,
            args.nda,
            args.target_summary_stats,]
    cmd.extend(args.training_summary_stats)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        sys.exit(stderr)
    reader = csv.DictReader(StringIO(stdout))
    result_rows = []
    for row in reader:
        row.update(extra_fields)
        result_rows.append(row)
    fieldnames = []
    fieldnames.extend([f for f in extra_fields.keys() if f not in fieldnames])
    fieldnames.extend([f for f in reader.fieldnames if f not in fieldnames])
    out = utility.open_output_file_for_csv_writer(
            filepath=args.output_filepath,
            append=args.append)
    writer = csv.DictWriter(out,
            fieldnames=fieldnames,
            lineterminator=os.linesep,
            )
    if not args.append:
        writer.writeheader()
    writer.writerows(result_rows)

if __name__ == "__main__":
    main()
