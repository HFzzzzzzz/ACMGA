import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-file", required=True, help="input file")
parser.add_argument("-o", "--output-file", required=True, help="output file")
parser.add_argument("-r", "--ref", required=True, help="value of ref")
parser.add_argument("-q", "--que", required=True, help="value of que")
args = parser.parse_args()
# with open(args.input_file) as input_file:

with open(args.input_file, "r", encoding='utf-8') as input_file, open(args.output_file, "w", encoding='utf-8') as output_file:
    for line in input_file:
        l = line.strip().split()  # noqa: E741
        lineone=l[0]
        linesix=l[5]

        # # Replace "chr10" at the beginning of the line with "id=foo|chr10".
        # Chr10 or chr10
        if('chr' in lineone):
            # chr10
            lineone_new = lineone.replace("chr", f"id={args.que}|chr", 1)
        elif ('Chr' in lineone):
            # Chr10
            lineone_new = lineone.replace("Chr", f"id={args.que}|Chr", 1)
        else:
            print("Error: please check the chromosome name.")
            break

        # Replace "chr10" in the sixth column with "id=bar|chr10".
        if('chr' in linesix):
            linesix_new = linesix.replace("chr", f"id={args.ref}|chr", 1)
        elif ('Chr' in linesix):
            linesix_new = linesix.replace("Chr", f"id={args.ref}|Chr", 1)
        else:
            print("Error: please check the chromosome name.")
            break
        l[0] = lineone_new
        l[5] = linesix_new
        line = "\t".join(l) + "\n"
        output_file.write(line)
