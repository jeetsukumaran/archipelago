#! /usr/bin/env python

import os
import collections
import argparse
import subprocess
import archipelago

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
    print(stdout.split("\n"))


if __name__ == "__main__":
    main()
