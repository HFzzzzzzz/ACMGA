import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-file", required=True, help="input file")
parser.add_argument("-o", "--output-file", required=True, help="output file")
parser.add_argument("-r", "--ref", required=True, help="value of ref")
parser.add_argument("-q", "--que", required=True, help="value of que")
args = parser.parse_args()
# with open(args.input_file) as input_file:

with open(args.input_file) as input_file, open(args.output_file, "w") as output_file:
    for line in input_file:
        lineone=line.split()[0]
        # print(lineone)
        linesix=line.split()[5]
        # print(linesix)
        # # Replace "chr10" at the beginning of the line with "id=foo|chr10".
        if('chr' in lineone):
            line = line.replace("chr", "id={}|chr".format(args.que), 1)
        # Replace "chr10" in the sixth column with "id=bar|chr10".
        if('chr' in linesix):
            line = line.replace("\tchr", "\tid={}|chr".format(args.ref), 1)
        if (line.split()[4] == "-"):
            start = int(line.split()[1]) - int(line.split()[3])
            end = start + int(line.split()[3]) - int(line.split()[2])
            start = str(start)
            end = str(end)
        else:
            start = line.split()[2]
            end = line.split()[3]
        spl = line.split()
        spl[2] = start
        spl[3] = end
        newline = '\t'.join(spl)
        output_file.write(newline)
        output_file.write("\n")
