import sys
import numpy
import pandas as pd
import pysam
import multiprocessing as mp
from collections import defaultdict
import argparse

def inputs_to_dict(input_queue, output_queue):
    while True:
        regions, myFileHandle = input_queue.get()
        if regions == "EMPTY":
            output_queue.put("DONE")
            break
        myDictionary = {}
        chr = regions[0]
        start = int(regions[1])
        end = int(regions[2])
        samfile = pysam.AlignmentFile(myFileHandle, "rb")
        for read in samfile.fetch(chr, start, end, multiple_iterators=True):            
            oldCigar = read.cigartuples
            if oldCigar == None:
                continue
            testTuple = tuple([read.reference_start, tuple(oldCigar), read.query_sequence])
            if testTuple in myDictionary:                
                myDictionary[testTuple] += 1        
            else:
                myDictionary[testTuple] = 1
        samfile.close()
        output_queue.put((regions, myFileHandle, myDictionary))

def filter_reads(myDictionary, start, end):
    start = start
    end = end
    length = end - start
    newDictionary = {}
    newList = []
    for a in myDictionary:
        pos, oldCigar, sequence = a
        thisRead = list(sequence)
        cigarPos = 0
        for cigarTuple in oldCigar:
            cigarType = int(cigarTuple[0])
            tupleSize = int(cigarTuple[1])
            if cigarType == 0 or cigarType == 7 or cigarType == 8: # alignment or sequence (mis)match
                cigarPos += tupleSize
            elif cigarType == 1: # insertion
                previous = thisRead[cigarPos-1]
                insertion = ''.join(thisRead[cigarPos:cigarPos+tupleSize])
                del thisRead[cigarPos:cigarPos+tupleSize]
                newString = str(previous) + "+" + str(tupleSize) + str(insertion)
                thisRead[cigarPos-1] = newString
            elif cigarType == 2: # deletion
                thisRead[cigarPos:cigarPos] = ["-"] * tupleSize
                cigarPos += tupleSize
            elif cigarType == 4: # soft-clipping
                del thisRead[cigarPos:cigarPos+tupleSize]
            elif cigarType == 3 or cigarType == 5 or cigarType == 6: # skipped region, hard-clipping, padding -> do nothing
                pass

        sequence = thisRead
        pos = int(pos)
        read_length = len(sequence)
        fin = pos + read_length
        lengthToAdd = 0
        new_start = 0
        new_fin = 0
        newRead = list()
        if pos <= start: 
            new_start = start - pos
            newRead = sequence[new_start:]
        if pos > start:
            lengthToAdd = pos - start
            newRead = [" "] * lengthToAdd + sequence
        if fin >= end:
            newRead = newRead[:length]
        if fin < end:
            lengthToAdd = end - fin
            newRead = newRead + [" "] * lengthToAdd
        totalCover = newRead.count(" ")
        if totalCover/float(length) > (1 - args_minCoverage):
            continue
        tupleToAppend = tuple(newRead)
        if tupleToAppend not in newDictionary:
            newDictionary[tupleToAppend] = myDictionary[a]
        else:
            newDictionary[tupleToAppend] += myDictionary[a]
        newList.append(newRead)
    return newDictionary, newList

def print_data(print_queue, final_queue):
    while True:
        toPrintSummary = []
        region, myDict = print_queue.get()
        if region == "EMPTY":
            final_queue.put("DONE")
            break
        chr = region[0]
        start = int(region[1]) - args_upstream
        end = int(region[2]) + args_downstream
        length = end - start
        toPrintSummary.append('{} {} {} {} {} {}'.format("Chromosome:", chr, "- Region:", start, "-", end))

        #convert alignment file to dictionary 
        controlDict = myDict[args_control]
        treatedDict = myDict[args_treated]

        filteredControlDict, filteredControlList = filter_reads(controlDict, start, end)
        filteredTreatedDict, filteredTreatedList = filter_reads(treatedDict, start, end)
        
        sampleDictList = [filteredControlDict, filteredTreatedDict]

        #get consensus sequence
        consensus = ""
        numpyData = numpy.array(filteredControlList)
        for column in numpyData.T:
            l = [x for x in column if x != " " and x != "-" and len(x) == 1]
            if not l:
                consensus += "N"
            else:
                (values,counts) = numpy.unique(l,return_counts=True)
                ind=numpy.argmax(counts)
                consensus += str(values[ind])

        if not consensus:
            consensus = "N" * length

        toPrintSummary.append('{} {}'.format("Consensus:", consensus))

        consensus = list(consensus)
        controlDict = {}
        mutDict = {}
        mutation_keys = set()

        total_sum_abs = 0
        isDeletion = False
        totalReadsList = []

        mutToRead = {}
        for count, d in enumerate(sampleDictList):
            mutationsDict = defaultdict(int)
            for sample in d:
                newMutation = []
                isDeletion = False
                myRead = list(sample)
                for i in range(len(sample)):
                    if sample[i] == ' ':
                        myRead[i] = consensus[i]
                    if sample[i] != consensus[i] and sample[i] != " ": #if mutation
                        if newMutation and newMutation[-1][0] == "-" and sample[i] == "-" and isDeletion: #if continuing a deletion
                            previous = newMutation[-1][1]
                            startIndex = str(previous.split(":")[0])
                            toAdd = startIndex + ":" + str(i)
                            newMutation[-1] = tuple(["-", toAdd])
                            isDeletion = True
                        else:
                            newMutation.append(tuple([sample[i], str(i)]))
                            if sample[i] == "-":
                                isDeletion = True
                            else:
                                isDeletion = False
                    else:
                        isDeletion = False
                newMutationTuple = tuple(newMutation)
                if newMutationTuple in mutationsDict:
                    mutationsDict[newMutationTuple] += d[sample]
                else:
                    mutationsDict[newMutationTuple] = d[sample]
                if newMutationTuple not in mutToRead:
                    mutToRead[newMutationTuple] = ''.join(myRead)
            total_reads = float(sum(mutationsDict.values()))
            totalReadsList.append(total_reads)
            if total_reads == 0:
                continue
            running_total = 0
            for a in sorted(mutationsDict.items(), key=lambda x: x[1], reverse = True):
                percentage_coverage = 100*a[1]/total_reads
                if a[0] not in controlDict:
                    controlDict[a[0]] = {}
                    controlDict[a[0]]['perc'] = 0
                    controlDict[a[0]]['count'] = 0
                if count == 0:
                    controlDict[a[0]]['perc'] = percentage_coverage
                    controlDict[a[0]]['count'] = a[1]
                if a[0] not in mutDict:
                    mutDict[a[0]] = {}
                    mutDict[a[0]]['perc'] = 0
                    mutDict[a[0]]['count'] = 0
                if count == 1:
                    mutDict[a[0]]['perc'] = percentage_coverage
                    mutDict[a[0]]['count'] = a[1]
            if count == 1:
                mutation_keys = set(mutDict.keys())
                mutation_keys.update(controlDict.keys())
                for a in mutation_keys:
                    total_sum_abs += abs(mutDict[a]['perc'] - controlDict[a]['perc'])
        toPrintSummary.append('{} {}'.format("Mutation Rate:", str(round(total_sum_abs/2, 3)) + "%"))
        toPrintSummary.append('{} {}'.format("Total Control Reads:", int(totalReadsList[0])))
        toPrintSummary.append('{} {}'.format("Total Treated Reads:", int(totalReadsList[1])))
        toPrintSummary.append("")
        toPrintSummary.append('{:35} {:30} {:>10} {:>10} {:>10} {:>10} {:>30}'.format("Sequence", "Mutation", "C-Cov", "C-#", "T-Cov", "T-#", "Difference"))
        toPrintData = []
        dataToPandas = []
        if mutation_keys:
            for a in mutation_keys:
                diff = round(mutDict[a]['perc']-controlDict[a]['perc'], 3)
                dataToPandas.append([mutToRead[a], ','.join([''.join(b) for b in a]), round(controlDict[a]['perc'], 3), controlDict[a]['count'], round(mutDict[a]['perc'], 3), mutDict[a]['count'], diff])
            pandasData = pd.DataFrame(dataToPandas, columns=["Sequence", "Mutations", "ContProp", "ContCount", "TreatProp", "TreatCount", "Difference"])
            pandasData = pandasData.values.tolist()
            formattedString = '{:35} {:30} {:10} {:10} {:10} {:10} {:30}'
            for a in pandasData:
                toPrintSummary.append(formattedString.format(*a))
        else:
            toPrintSummary.append("No mutations")
        toPrintSummary.append("")
        separator = "-" * 50    
        toPrintSummary.append(separator)
        toPrintSummary.append("")
        final_queue.put(toPrintSummary)


parser = argparse.ArgumentParser()
parser.add_argument('regions', help="List of regions (bed)")
parser.add_argument('control', help="Control BAM file")
parser.add_argument('treated', help="Treated BAM file")
parser.add_argument('-up', type=int, nargs='?', help="Extend all regions upstream by X positions", default=0, const=0)
parser.add_argument('-down', type=int, nargs='?', help="Extend all regions downstream by X positions", default=0, const=0)
parser.add_argument('-mc', type=float, nargs='?', help="Minimum read coverage over region (decimal between 0 and 1)", default=0.5, const=0.50)
parser.add_argument('-mr', type=float, nargs='?', help="Minimum relative read frequency to classify as significant mutation (percentage between 0 and 100)", default=1, const=1)
parser.add_argument('-o', nargs='?', help="Optional Output File")

args = parser.parse_args()
args_regions = args.regions
args_control = args.control
args_treated = args.treated
args_upstream = args.up
args_downstream = args.down
args_minCoverage = args.mc
args_minFreq = args.mr
args_out = args.o

regions = []
with open(args_regions) as f:
    for line in f:
        words = line.split()
        chrom, start, stop = words[0:3]
        if chrom[0:3] == "chr":
            chrom = chrom[3:]
        regions.append(tuple([chrom, start, stop]))

print("Parsing Input")

input_queue = mp.Queue()
output_queue = mp.Queue()
print_queue = mp.Queue()
final_queue = mp.Queue()
print_processes = []
num_threads = mp.cpu_count()
processes = [mp.Process(target=inputs_to_dict, args=(input_queue, output_queue)) for x in range(num_threads)]
for p in processes:
    p.start()
num_done = 0

parsedData = {}
for r in regions:
    parsedData[r] = {}

for region in regions:
    input_queue.put((region, args_control))
    input_queue.put((region, args_treated))
for x in range(num_threads):
    input_queue.put(("EMPTY", "EMPTY"))

print("Analysing Data")
while num_done < num_threads:
    output = output_queue.get()
    if output == "DONE":
        num_done += 1
        p = mp.Process(target=print_data, args=(print_queue, final_queue))
        print_processes.append(p)
        p.start()
    else:
        region = output[0]
        myFile = output[1]
        myDict = output[2]
        parsedData[region][myFile] = myDict
        if len(parsedData[region]) == 2: #both treated and control have been added
            print_queue.put((region, parsedData[region]))
for p in processes:
    p.join()

for x in range(num_threads):
    print_queue.put(("EMPTY", "EMPTY"))

num_done = 0

dataToPrint = {}
while num_done < num_threads:
    output = final_queue.get()
    if output == "DONE":
        num_done += 1
    else:
        dataToPrint[output[0]] = '\n'.join(output)

if args_out:
    f = open(args_out, 'w')
    for key, value in sorted(dataToPrint.items(), key=lambda x: x[0]): 
        f.write(value)
    f.close()
else:
    for key, value in sorted(dataToPrint.items(), key=lambda x: x[0]): 
        print(value)

for p in print_processes:
    p.join()