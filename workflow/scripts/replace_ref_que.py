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
        lineone=line.split()[0]
        # print(lineone)
        linesix=line.split()[5]
        # print(linesix)
        # # Replace "chr10" at the beginning of the line with "id=foo|chr10".
        # Chr10 or chr10
        if('chr' in lineone):
            # chr10
            line = line.replace("chr", "id={}|chr".format(args.que), 1)
        elif ('Chr' in lineone):
            # Chr10
            line = line.replace("Chr", "id={}|Chr".format(args.que), 1)
        else:
            print("Error: please check the chromosome name.")
            break

        # Replace "chr10" in the sixth column with "id=bar|chr10".
        if('chr' in linesix):
            line = line.replace("\tchr", "\tid={}|chr".format(args.ref), 1)
        elif ('Chr' in linesix):
            line = line.replace("\tChr", "\tid={}|Chr".format(args.ref), 1)
        else:
            print("Error: please check the chromosome name.")
            break
        output_file.write(line)
