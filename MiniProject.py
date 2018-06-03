#!/usr/bin/env python
"""Simple script to visualize the coverage of some gene regions.

Requires a gene_models file (i.e., gene_models.txt) and a data_file
(i.e., input_data.txt). If anyone actually wants to try the script out, let me know and
I'll pass along the files.

By default draws two columns of nine rows each, generating a single plot for
each gene.
"""

import argparse
import matplotlib.pyplot as plt


def arg_parser():
    """Parse the arguments passed when calling the script."""
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene_models", required = True,
                        help = "Gene Models file")
    parser.add_argument("-d", "--data", required = True, help = "Data file")
    parser.add_argument("--not_all_genes", required = False,
                        action = "store_true",
                        help = "Only display genes for which a sample has data")
    parser.add_argument("--draw_exons", required = False, action = "store_true",
                        help = "Draw gene regions relative to the exons")
    parser.add_argument("--write_strand", required = False, action = "store_true",
                        help = "Writes the genomic and read strands")
    args = parser.parse_args()
    return args


def read_models_file():
    """Read the gene_models file and extracts the data."""
    global gene_list  # List of all the genes
    global gene_dict  # Dictionary of all genes and their exons
    gene_list = []
    gene_dict = {}
    models_file = open(args.gene_models, "r")
    next(models_file)  # Skip the first (header) line of the models_file

    # Reads through the models_file and pulls out each gene and its exon coordinates
    for line in models_file:
        line = line.strip().split("\t")
        gene = line[0]
        exon_start = int(line[4])
        exon_end = int(line[5])
        exon_num = str(line[6])
        if gene not in gene_list:
            gene_list.append(gene)
            gene_dict[gene] = {}
        gene_dict[gene][exon_num] = (exon_start, exon_end)
    models_file.close()
    return(gene_list, gene_dict)


def read_data_file():
    """Read the data_file and extract the information."""
    global sample_list  # List of all samples
    global sample_dict  # Dictionary of all samples with their genes and gene regions
    global gene_region_dict  # Dictionary of all genes and their regions
    sample_list = []
    sample_dict = {}
    gene_region_dict = {}
    data_file = open(args.data, "r")
    next(data_file)  # Skip the first (header) line

    # Reads through the data_file and pulls out samples, genes, gene_regions,
    # coordinates, and strands
    for line in data_file:
        line = line.strip().split("\t")
        sample = line[0]
        gene_region = line[1]
        gene = line[2]
        tx_strand = line[6]
        read_strand = line[7]
        start = int(line[8])
        end = int(line[9])
        templates = int(line[10])

        # If the sample isn't in the sample_list, add it
        if sample not in sample_list:
            sample_list.append(sample)
            sample_dict[sample] = {}

        # If the gene isn't in the sample, add it
        if gene not in sample_dict[sample]:
            sample_dict[sample][gene] = {}
        sample_dict[sample][gene][gene_region] = [tx_strand, read_strand,
                                                  start, end, templates]

        # If the gene isn't in the gene_region_dict, add it
        if gene not in gene_region_dict:
            gene_region_dict[gene] = []

        # If the gene_region isn't in the gene, add it
        if gene_region not in gene_region_dict[gene]:
            gene_region_dict[gene].append(gene_region)
    data_file.close()
    return(sample_list, sample_dict, gene_region_dict)


def exon_coords(gene):
    """Calculate the max and min exon coordinates for a given gene."""
    # exon_min and exon_max set arbitrary Min and Max exon coordinates
    exon_min = gene_dict[gene]["1"][0]
    exon_max = gene_dict[gene]["1"][0]

    # Determines the actual Min and Max exon coordinates for the gene
    for exon in gene_dict[gene]:
        for coord in gene_dict[gene][exon]:
            if coord < exon_min:
                exon_min = coord
            elif coord > exon_max:
                exon_max = coord
    return(exon_min, exon_max)


def draw_gene_regions(sample, gene, gene_plot):
    """Draw the sequenced regions.

    Draws each sequenced region for the gene. Height of the bars correspond to the
    number of Templates.

    If the number of templates is 0 < templates < 200, draws the bars orange,
    otherwise they are grey.

    Also labels each sequenced region, and draws (+/-) for the strands if the
    --write_strand option is turned on.
    """
    num_regions = len(gene_region_dict[gene])
    regions_scaler = 1000 / num_regions  # Scales the width of each region so they fit in the subplot
    for x in range(num_regions):
        region = gene_region_dict[gene][x]

        # Only draws bars if the sample has data for the region, otherwise
        # draws a line (really, a bar with 0.1 height)
        if region in sample_dict[sample][gene]:
            coverage = sample_dict[sample][gene][region][4]
            if coverage >= 200:
                gene_plot.bar(regions_scaler * x, coverage, regions_scaler, 1,
                              align = "edge", linewidth = 1, log = True,
                              color = "#BBB4AD")
            # If the region has fewer than 200 temples, draws the box orange
            # and prints a warning
            else:
                gene_plot.bar(regions_scaler * x, coverage, regions_scaler, 1,
                              align = "edge", linewidth = 1, log = True,
                              color = "#FF7A00")
                print("Warning: {} in sample {} had fewer than 200 templates"
                      .format(region, sample))
        else:
            gene_plot.bar(regions_scaler * x, 0.1, regions_scaler, 1,
                          align = "edge", linewidth = 1, log = True,
                          color = "#BBB4AD")

        # Writes the region name (everything after the last '_')
        text_x = (2 * (regions_scaler * x) + regions_scaler) / 2
        plt.text(text_x, 1.1, region[region.rfind("_") + 1:], fontsize = 6,
                 rotation = 90, ha = "center", va = "bottom")

        # If '--write_strand' is turned on, writes (+/-) for the transcript and
        # read orientation relative to the genomic reference
        if args.write_strand is True and region in sample_dict[sample][gene]:
            strand_text = "{}{}".format(sample_dict[sample][gene][region][0],
                                        sample_dict[sample][gene][region][1])
            plt.text(text_x, coverage + (coverage * 0.1), strand_text,
                     fontsize = 6, ha = "center", va = "bottom")


def draw_exons(gene, gene_plot):
    """Draw exons for the gene."""
    exon_min, exon_max = exon_coords(gene)
    exon_scaler = exon_max - exon_min  # Determines the total length of all exons
    nt_length = 1000 / exon_scaler  # Determines the width of each nucleotide, so they fit in the subplot

    # Draws each exon in blue
    for exon in gene_dict[gene]:
        coords = [nt_length * (exon_max - x) for x in gene_dict[gene][exon]]
        gene_plot.plot(coords, [1.1, 1.1], linewidth = 1, color = "#07B8FB")


def draw_coverage(sample, gene, gene_plot):
    """Draw coverage for the sequenced region, relative to its position by exons."""
    exon_min, exon_max = exon_coords(gene)
    exon_scaler = exon_max - exon_min  # Determines the total length of all exons
    nt_length = 1000 / exon_scaler  # Determines the width of each nucleotide, so they fit in the subplot

    # Draws a bar for every region in the gene, with its height based on the
    # sequencing coverage
    for region in sample_dict[sample][gene]:
        template_coords = sample_dict[sample][gene][region][2:4]
        template_coords = [nt_length * (exon_max - x) for x in template_coords]
        coverage = sample_dict[sample][gene][region][4]
        if coverage >= 200:
            gene_plot.bar(template_coords[0], coverage, abs(template_coords[0] - template_coords[1]),
                          1, align = "edge", linewidth = 0, log = True,
                          color = "#BBB4AD", zorder = 0)

        # If the region has fewer than 200 temples, draws the box orange
        # and prints a warning
        else:
            gene_plot.bar(template_coords[0], coverage, abs(template_coords[0] - template_coords[1]),
                          1, align = "edge", linewidth = 0, log = True,
                          color = "#FF7A00", zorder = 10)
            print("Warning: {} in sample {} had fewer than 200 templates"
                  .format(region, sample))


def start_drawing(sample, num_genes, genes_to_draw):
    """Format the figure and call the drawing functions for each gene."""
    gene_count = 1
    # For each gene, makes a subplot, formats it, and calls the appropriate
    # draw function(s)
    for gene in genes_to_draw:
        gene_plot = plt.subplot((num_genes / 2) + (num_genes % 2), 2, gene_count)
        plt.yscale("log")
        plt.axis([0, 1000, 1, 100000])
        plt.xlabel(gene)
        plt.ylabel("Templates", fontsize = 8)
        gene_plot.spines['top'].set_visible(False)
        gene_plot.spines['bottom'].set_visible(False)
        gene_plot.spines['right'].set_visible(False)
        gene_plot.tick_params(labelbottom = False, bottom = False)

        # If the gene is present in the sample, draw it. Otherwise just draw
        # a line
        if gene in sample_dict[sample]:
            if args.draw_exons is True:
                draw_exons(gene, gene_plot)
                draw_coverage(sample, gene, gene_plot)
            else:
                draw_gene_regions(sample, gene, gene_plot)
        else:
            gene_plot.bar(0, 0.1, 1000, 1,
                          align = "edge", linewidth = 1, log = True,
                          color = "#BBB4AD")

        plt.yticks([10, 100, 1000, 10000, 100000], fontsize=8)
        gene_count += 1


def run_script():
    """Call all the functions of the script."""
    args = arg_parser()
    gene_list, gene_dict = read_models_file()
    sample_list, sample_dict, gene_region_dict = read_data_file()

    # Makes a figure for each sample in the data_file
    for sample in sample_list:
        print("############################################################")
        print("Working on {}".format(sample))
        num_genes = len(gene_list)
        genes_to_draw = gene_list

        # If not_all_genes is turned on, sets num_genes and genes_to_draw to
        # only include genes present in the sample
        if args.not_all_genes is True:
            num_genes = len(sample_dict[sample])
            genes_to_draw = sample_dict[sample]
        plt.figure(figsize = (8.5, 11))
        start_drawing(sample, num_genes, genes_to_draw)  # Makes the magic happen
        plt.tight_layout()
        print("Saving figure as {}.pdf\n\n".format(sample))
        plt.savefig("{}.pdf".format(sample))
        plt.clf()  # Clears the figure so bars don't carry over between samples


if __name__ == "__main__":
    run_script()
    print("All done!")
