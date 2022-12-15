import argparse
import sys
import os

def localize(args):
    import localize/localize_pseudogenes
    import localize/find_unique_sites

    pseudogene_dict = localize_pseudogenes.create_pseudogene_dict(args.summary_file)
    gene_dict = localize_pseudogenes.create_gene_info_dict(args.gene_pos)

    localize_pseudogenes.pseudogene_localization(args.VCF, pseudogene_dict, gene_dict, args.out_path, args.prefix, args.r)

    find_unique_sites.unique_pos(args.out_path, args.Unique_len, args.prefix)

def stat(args):
    import statistics/binomial_analysis

    MC_list = binomial_analysis.insert_generation(args.n, args.output_file, args.out_path, args.prefix)
    Actual_list = binomial_analysis.loc_file_to_list(args.insert_locations)

    binomial_analysis.chromosome_binomial(MC_list,Actual_list, args.prefix, args.out_path)
    binomial_analysis.check_all_track(args.track_file, Actual_list, MC_list, args.output_file, args.out_path, args.prefix)


def main():
    parser = argparse.ArgumentParser(description= "Discover and localize non-reference Pseudogenes and run statistics on them.")
    subparsers = parser.add_subparsers()

    parser_localize = subparsers.add_parser("localize", help = "localize help")
    parser_localize.add_argument("summary_file", type=str, help = "Path to the summary file created in the Discover step of the Processen pipeline.")
    parser_localize.add_argument("gene_pos", type = str, help = "Path to a bed file containing gene positioning information in a tab separated format. "
                                                                    "Each line should be formatted as <gene> <start> <end> <chromosome>")
    parser_localize.add_argument("VCF", type = str, help = "Path to a file containing the paths to the VCF file that has sequencing information of the individual of interest.")
    parser_localize.add_argument("out_path", type = str, help = "Path to the folder where output files (files contianing pseudogene positioning information for each individual) should be created. The directory should be empty. Should end on /.")
    parser_localize.add_argument("prefix", type = str, help = "Prefix for the output files. The final output file name will be <prefix>_unique_inserts.txt, it will be created in the current folder.")
    parser_localize.add_argument("--Unique_len", "-l", type = int, default = 100, help = "The minimum distance between two insert sites to be considered unique")
    parser_localize.add_argument("--r", type=int, default=1000, help = "the max distance a point can be from the o.g. location of the gene to be considered further. If this number is made smaller less variant calls will be attributed to this gene. If it is made too big, variants that are outside of this gene may also be considered for the location of the pseudogene.")
    parser_localize.set_defaults(func=localize)

    parser_statistics = subparsers.add_parser("statistics", help = "statistics help")
    parser_statistics.add_argument("prefix", type=str, help = "Prefix for the output files.")
    parser_statistics.add_argument("insert_locations", type = str, help = "Path to the file with unique insert sites generated in the localize step of the Processen pipeline.")
    parser_statistics.add_argument("track_file", type = str, help = "Path to a file containing the path to a bed file containing locations of each track that should be analysed for overrepresentation."
    parser_statistics.add_argument("--out_path", type = str, default = '', "Path to directory where the output files should be stored. Should end on /.")
    parser_statistics.add_argument("--n", type = int, default = 100000, help = "Number of random insert sites to generate for binomial analysis.")
    parser_statistics.add_argument("--output_file", type = bool, default = False, help = "If output files for all intermediate steps should be made put True, else put False.")
    parser_localize.set_defaults(func = stat)

    args = parser.parse_args(sys.argv[1:])
    args.func(args)

if __name__ == '__main__':
    main()