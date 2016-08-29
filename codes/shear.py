 # -*- coding: utf-8 -*-
"""
SHEAR: Simple Handler for Error and Adapter Removal
James B. Pease
Version Alpha 006

VERSION HISTORY:

v.003 Alpha Release
v.004 Added support for combining multiple pairs of input files
v.005 Major update, fixed filtered output, general clean-up, removed GC filter
v.006 Minor fixes to default parameters
"""
SHEARVERSION = "Alpha 005"

import sys, argparse, os, subprocess
from time import time
from datetime import datetime
from math import log

STAT_STRING = ("""
original reads:\t{!s}
original bases:\t{!s}
reads with adapters removed:\t{!s}
reads partially trimmed for quality:\t{!s}\t{!s}
bases trimmed for quality:\t{!s}\t{!s}
reads filtered for length:\t{!s}\t{!s}
bases filtered for length:\t{!s}\t{!s}
reads filtered for low avg quality:\t{!s}\t{!s}
bases filtered for low avg quality:\t{!s}\t{!s}
reads filtered for ambiguity:\t{!s}\t{!s}
bases filtered for ambiguity:\t{!s}\t{!s}
reads filtered for low info:\t{!s}\t{!s}
bases filtered for low info:\t{!s}\t{!s}
reads filtered because unpaired:\t{!s}\t{!s}
bases filtered because unpaired:\t{!s}\t{!s}
reads removed:\t{!s}\t{!s}
reads trimmed:\t{!s}\t{!s}
bases removed:\t{!s}\t{!s}""")


class TrimStat(object):
    """Trimming Statistics Handler"""

    def __init__(self):
        # Quartets represent [FQ1reads, FQ1bases, FQ2reads, FQ2bases]
        self.stats = dict(tuple([x, [0, 0, 0, 0]]) for x in [
            'scythe', 'total', 'ambig', 'totalqual',
            'trim', 'lowinfo', 'unpaired', 'tooshort'])

    def p1add(self, statname, value):
        """Add values to p1 stats"""
        self.stats[statname][0] += 1
        self.stats[statname][1] += value
        return ''

    def p2add(self, statname, value):
        """Add values tp p2 stats"""
        self.stats[statname][2] += 1
        self.stats[statname][3] += value
        return ''

    def p1scythe(self, value):
        """Add to scythe stat for p1"""
        self.stats['scythe'][0] += value

    def p2scythe(self, value):
        """Add to scythe stat for p2"""
        self.stats['scythe'][2] += value

    def stat_list(self, pair1=False, pair2=False, ndigits=3):
        """Create statistics output strings"""
        reads_rem = 0
        bases_rem = 0
        reads_trim = 0
        (i, j) = pair1 and (0, 1) or pair2 and (2, 3)
        xlist = [self.stats['total'][i],
                 self.stats['total'][j],
                 self.stats['scythe'][i]]

        for statname in ['trim', 'tooshort', 'totalqual', 'ambig',
                         'lowinfo', 'unpaired']:
            xlist.extend([self.stats[statname][i],
                          zdiv(self.stats[statname][i],
                               self.stats['total'][i], ndigits),
                          self.stats[statname][j],
                          zdiv(self.stats[statname][j],
                               self.stats['total'][j], ndigits)])
            if statname in ['trim', 'scythe']:
                reads_trim += self.stats[statname][i]
            else:
                reads_rem += self.stats[statname][i]
            bases_rem += self.stats[statname][j]
        xlist.extend([reads_rem, zdiv(reads_rem,
                                      self.stats['total'][i], ndigits),
                      reads_trim, zdiv(reads_trim,
                                       self.stats['total'][i], ndigits),
                      bases_rem, zdiv(bases_rem,
                                      self.stats['total'][j], ndigits)])
        return STAT_STRING.format(*xlist)

def timestamper():
    """Returns dash-separated timestamp string"""
    return datetime.now().strftime('%Y_%m_%d_%H_%M_%S_%f')

def zdiv(num, den, ndigits=-1):
    """Calculate a fraction, but avoid zero denominator errors"""
    if not den:
        return 0.0
    else:
        num = float(num) / float(den)
        if ndigits != -1:
            num = round(num, ndigits)
        return num

def complement(seq):
    """Make complement of DNA sequence"""
    seq = seq[::-1].lower().replace('a', 'T')
    return seq.replace('t', 'A').replace('g', 'C').replace('c', 'G')

def mutual_info(seq):
    """Calculate a simple dinucleotide mutual information score"""
    if len(seq) < 3:
        return 0.0
    dinucs = {}
    monucs = {}
    for j in xrange(len(seq) - 1):
        dinucs[seq[j:j+2]] = dinucs.get(seq[j:j+2], 0) + 1
        monucs[seq[j]] = monucs.get(seq[j], 0) + 1
    monucs[seq[-1]] = monucs.get(seq[-1], 0) + 1
    total_monucs = float(sum(monucs.values()))
    total_dinucs = float(sum(dinucs.values()))
    mutual_information = 0
    if not total_dinucs:
        return 0.0
    for dinuc in dinucs:
        p_dinuc = float(dinucs[dinuc]) / total_dinucs
        mutual_information += (p_dinuc * log(p_dinuc / (
            (float(monucs[dinuc[0]]) / total_monucs) *
            (float(monucs[dinuc[1]]) / total_monucs))))
    return mutual_information

def run_scythe(fastqpath, adapters=None, out=None, **kwargs):
    """Runs the Scythe Bayesian Trimmer"""
    scythelog = open(kwargs.get('logfilepath', 'scythe.log'), "w")
    quality = kwargs.get('quality', 'sanger')
    if quality in ['illumina_1.3', 'illumina_1.5']:
        quality = 'illumina'
    elif quality in ['illumina_1.8']:
        quality = 'sanger'
    arr_cmd = [kwargs.get('execpath', 'scythe'),
               '-a', adapters, '-q', quality,
               '-p', str(kwargs.get('prior', 0.1)),
               '-o', out,
               '-n', str(kwargs.get('min_match', 5)),
               '-M', str(kwargs.get('min_keep', 0)),
               fastqpath]
    if kwargs.get('verbose', True):
        print("executing subprocess {!s}").format(' '.join(arr_cmd))
    process = subprocess.Popen(' '.join(arr_cmd),
                               shell=True, stderr=scythelog, stdout=scythelog)
    process.communicate()
    contaminated = 0
    scythelog = open(kwargs.get('logfilepath', 'scythe.log'), 'r')
    line = scythelog.readline()
    while line:
        if line.startswith('contaminated'):
            contaminated = int(line[line.find(': ') + 2:line.find(', ')])
            break
        line = scythelog.readline()
    scythelog.close()
    return contaminated

def detect_fastq_format(filepath, ntest=1000000):
    """Determine FASTQ Format"""
    i = 0
    with open(filepath) as fastqfile:
        code = ''
        while code in ['', 'phred33', 'phred64'] and i < ntest:
            i += 1
            fastqfile.readline() #Skip Header
            bases = fastqfile.readline()
            fastqfile.readline() #Skip spacer
            quals = fastqfile.readline()
            bases = bases.rstrip()
            if code == 'phred33':
                if 'J' in quals:
                    code = "illumina_1.8"
                elif any(x in quals for x in '!"'):
                    code = "sanger"
            elif code == 'phred64':
                if any(x in quals for x in ';<=>?'):
                    code = "solexa"
                elif any(x in quals for x in '@"'):
                    code = 'illumina13'
            else:
                if any(x in quals for x in '!"#$%&\'()*+,-./0123456789:'):
                    code = 'phred33'
                    if 'J' in quals:
                        code = "illumina_1.8"
                    elif any(x in quals for x in '!"'):
                        code = "sanger"
                elif any(x in quals for x in 'KLMNOPQRSTUVWXYZ[]^_`abcdefgh'):
                    code = 'phred64'
                    if any(x in quals for x in ';<=>?'):
                        code = "solexa"
                    elif any(x in quals for x in '@"'):
                        code = 'illumina13'
        while i < ntest:
            fastqfile.readline()
            bases = fastqfile.readline()
            fastqfile.readline()
            fastqfile.readline()
            if not bases:
                break
            bases = bases.rstrip()
    if code == 'phred64':
        code = "illumina_1.5"
    elif code == 'phread33':
        code = "sanger"
    return code

def generate_adapters(outputpath, mode='', fwd_barcodes=None,
                      rev_barcodes=None):
    """Generate Adapter File for use with Scythe"""
    if not fwd_barcodes:
        fwd_barcodes = ['']
    if not rev_barcodes:
        rev_barcodes = ['']
    if mode in ['TruSeqDualIndex']:
        fwd_adapt = [("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[barcode]"
                      "ATCTCGTATGCCGTCTTCTGCTTGAAAAA"),
                     ("ATATCGGAAGAGCACACGTCTGAACTCCAGTCAC[barcode]"
                      "ATCTCGTATGCCGTCTTCTGCTTGAAAAA"),
                     ("TTTTTCAAGCAGAAGACGGCATACGAGAT[barcode]"
                      "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")]
        rev_adapt = [("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT[barcode]"
                      "GTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA"),
                     ("ATATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT[barcode]"
                      "GTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA"),
                     ("TTTTTAATGATACGGCGACCACCGAGATCTACAC[barcode]"
                      "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")]
    elif mode in ['TruSeq']:
        fwd_adapt = [("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[barcode]"
                      "ATCTCGTATGCCGTCTTCTGCTTGAAAAA"),
                     ("ATATCGGAAGAGCACACGTCTGAACTCCAGTCAC[barcode]"
                      "ATCTCGTATGCCGTCTTCTGCTTGAAAAA"),
                     ("TTTTTCAAGCAGAAGACGGCATACGAGAT[barcode]"
                      "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")]
        rev_adapt = [("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"
                      "ATCTCGGTGGTCGCCGTATCATTAAAAA"),
                     ("ATATCGGAAGAGCGTCGTGTAGGGAAAGAG"
                      "TGTGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA"),
                     ("TTTTTAATGATACGGCGACCACCGAGAT"
                      "CTACACACACTCTTTCCCTACACGACGCTCTTCCGATCT")]
    with open(outputpath, 'w') as outfile:
        for k, fwd_barcode in enumerate(fwd_barcodes):
            fwd_comp = complement(fwd_barcode)
            outfile.write((">TruSeq_fwd_adapter_{!s}\n{!s}\n"
                           ">TruSeq_fwd_adapter_alt_{!s}\n{!s}\n"
                           ">TruSeq_fwd_adapter_comp_{!s}\n{!s}\n"
                          ).format(k+1, fwd_adapt[0].replace('[barcode]',
                                                             fwd_barcode),
                                   k+1, fwd_adapt[1].replace('[barcode]',
                                                             fwd_barcode),
                                   k+1, fwd_adapt[2].replace('[barcode]',
                                                             fwd_comp)))
        for k, rev_barcode in enumerate(rev_barcodes):
            rev_comp = complement(rev_barcode)
            outfile.write((">TruSeq_rev_adapter_{!s}\n{!s}\n"
                           ">TruSeq_rev_adapter_alt_{!s}\n{!s}\n"
                           ">TruSeq_rev_adapter_comp_{!s}\n{!s}\n"
                          ).format(k+1, rev_adapt[0].replace('[barcode]',
                                                             rev_comp),
                                   k+1, rev_adapt[1].replace('[barcode]',
                                                             rev_comp),
                                   k+1, rev_adapt[2].replace('[barcode]',
                                                             rev_barcode)))
    return ''



def main(arguments=sys.argv[1:]):
    """FASTQ Trimmer Main Function"""
    time0 = time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq1", required=True, nargs='*')
    parser.add_argument("--fq2", nargs='*')
    parser.add_argument("--out1", required=True)
    parser.add_argument("--out2")
    parser.add_argument("--platform", choices=['TruSeq', 'TruSeqDualIndex'],
                        default='TruSeq', help="Sequencing Platform")
    parser.add_argument("--barcodes1", nargs='*',
                        help="One or more barcodes for the first reads")
    parser.add_argument("--barcodes2", nargs='*',
                        help="One or more barcodes for the second reads" +
                        "paired-end TruSeqDualIndex only")
    parser.add_argument("--trimfixed", default="0:0",
                        help="remove fixed number of bases from the FRONT:END")
    parser.add_argument("--trimqual", default="20:20",
                        help=("remove bases below this quality score the"
                              + " FRONT:END"))
    parser.add_argument("--trimqualpad", default="0:0",
                        help=("removing this many extra bases next to a "
                              "low-quality base from FRONT:END"))
    parser.add_argument("--trimpolyat", type=int, default=12,
                        help=("remove a poly-A or poly-T repeat of"
                              "at least this length from the front or end"))
    parser.add_argument("--trimambig", action="store_true", default=True,
                        help="remove N's from both ends")
    parser.add_argument("--filterlength", type=int, default=1,
                        help="minimum read length")
    parser.add_argument("--filterquality", type=int, default=3,
                        help=("filter reads with pre-trimming mean quality"
                              "below this score"))
    parser.add_argument("--filterunpaired", action="store_true",
                        help="remove unpaired reads (paired mode only)")
    parser.add_argument("--filterlowinfo", type=float, default=0.,
                        help=("filter out reads with low information"
                              "(default=0.5 dinucleotide MI score,"
                              "set to 0 to skip)"))
    parser.add_argument("--filterambig", type=int, default=5,
                        help=("filter reads with more than this number"
                              "of ambiguous bases, (default=2, set=0 to skip"))
    parser.add_argument("--trimpattern3", nargs='*',
                        help="remove these specific base sequences from end")
    parser.add_argument("--skipscythe", action="store_true",
                        help="skip scythe 3' adapter removal")
    parser.add_argument("--execscythe",
                        help="path of the scythe executable")
    parser.add_argument("--scytheprior", type=float, default=0.1)
    parser.add_argument("--log", help="specify log file path"
                        + "default is shear_TIMESTAMP")
    parser.add_argument("--cleanup", choices=["none", "tempfastq",
                                              "exceptadapters", "all"],
                        default="all",
                        help=("choose which temporary files to remove ("
                              "none=keep all; "
                              "tempfastq=remove temporary fastq files"
                              "from scythe; "
                              "exceptadapters=remove temp fastq and log"
                              "files from skype but keep adapter file; "
                              "all=keep all temporary files;)"))
    parser.add_argument("--filterfile",
                        help="specify output file prefix for filtered reads")
    parser.add_argument("--tempdir",
                        help="directory to use for temporary files")
    parser.add_argument("--cleanheader", action="store_true",
                        help=("removes any additional terms from the"
                             "header line (useful after STAR)"))
    args = parser.parse_args(args=arguments)

#### ====== Establish Log and Temporary Directory ===========
    timestamp = timestamper()
    if not args.tempdir:
        args.tempdir = '.'
    args.tempdir = os.path.abspath(args.tempdir)
    if not os.path.exists(args.tempdir):
        os.mkdir(args.tempdir)
    if not args.log:
        args.log = "shear_" + timestamp + ".log"
    args.log = os.path.abspath(args.log)
    mainlog = open(args.log, 'w')
    mainlog.write("shear version={!s}\n{!s}\n".format(SHEARVERSION,
                                                      " ".join(sys.argv)))

#### ====== Additional Argument Checks ============
    try:
        trimfixed = tuple([int(x) for x in args.trimfixed.split(':')][0:2])
    except:
        msg = "--trimfixed must be in form FRONT:END"
        raise SyntaxError(msg)
    try:
        trimqual = tuple([int(x) for x in args.trimqual.split(':')][0:2])
    except:
        msg = "--trimqual must be in form FRONT:END"
        mainlog.write(msg + "\n")
        raise SyntaxError(msg)
    try:
        trimqualpad = tuple([int(x) for x in args.trimqualpad.split(':')][0:2])
    except:
        msg = "--trimqualpad must be in form FRONT:END"
        mainlog.write(msg + "\n")
        raise SyntaxError(msg)
    polynucs = []
    if args.trimambig:
        polynucs.append(('N',1))
    if args.trimpolyat:
        polynucs.extend([('A', args.trimpolyat), ('T',args.trimpolyat)])

#### ====== Establish file paths and paired/single mode =======
    args.fq1 = [os.path.abspath(x) for x in args.fq1]
    args.out1 = os.path.abspath(args.out1)
    outfq1 = open(args.out1, 'w')
    if args.filterfile:
        filterfq1 = open(os.path.abspath("{!s}_p1.fastq".format(
            args.filterfile)), 'w')
    paired_end = False
    if args.fq2:
        args.fq2 = [os.path.abspath(x) for x in args.fq2]
        paired_end = True
        if not args.out2:
            msg = "--fq2 specified without --out2"
            mainlog.write(msg + "\n")
            raise RuntimeError(msg)
        args.out2 = os.path.abspath(args.out2)
        if args.filterfile:
            filterfq2 = open(os.path.abspath("{!s}_p2.fastq".format(
                args.filterfile)), 'w')
        outfq2 = open(args.out2, 'w')
        input_fq = zip(args.fq1, args.fq2)
    else:
        input_fq = zip(args.fq1, ['']*len(args.fq1))
#### ===== Establish stats container ==========
    stats = TrimStat()
#### ==== Scythe Adapter Preparation =======
    if args.skipscythe:
        mainlog.write("Skipping Scythe...\n")
    else:
        mainlog.write("Using temporary working directory: {!s}\n".format(
            args.tempdir))
        adapterpath = "{!s}/scythe_adapters_{!s}.fasta".format(args.tempdir,
                                                               timestamp)
        generate_adapters(adapterpath, mode=args.platform,
                          fwd_barcodes=args.barcodes1,
                          rev_barcodes=args.barcodes2)
        mainlog.write("Adapters written to: {!s}".format(adapterpath))

#### ======= BEGIN ITERATION =============
    ndex = 0
    tempfile1, tempfile2 = "", ""
    for (fq1, fq2) in input_fq:
        code_p1 = detect_fastq_format(fq1,)
        if not code_p1:
            msg = "fq1 quality score scale cannot be determined"
            mainlog.write(msg + "\n")
            raise RuntimeError(msg)
        else:
            mainlog.write("fq1 quality scale detected: {!s}\n".format(code_p1))
        if paired_end:
            code_p2 = detect_fastq_format(fq2)
            if not code_p2:
                msg = "fq2 quality score scale cannot be determined"
                mainlog.write(msg + "\n")
                raise RuntimeError(msg)
            else:
                mainlog.write("fq2 quality scale detected: {!s}\n".format(
                    code_p2))
            if code_p1 != code_p2:
                msg = ("fq1 and fq2 quality score scales do not match "
                       "({!s}, {!s})").format(code_p1, code_p2)
                mainlog.write(msg + "\n")
                raise RuntimeError(msg)
######## ==== Scythe Adapter Trimming on Pair 1 =======
        if not args.skipscythe:
            tempfile1 = os.path.abspath("{!s}/scythetemp_p1_{!s}_{!s}"
                                        ".fastq".format(
                                            args.tempdir, timestamp, ndex))
            logfile1 = os.path.abspath("{!s}/scythe_p1_{!s}_{!s}.log".format(
                args.tempdir, timestamp, ndex))
            mainlog.write("Running Scythe on fq1...\n")
            stats.p1scythe(run_scythe(fq1, adapters=adapterpath, out=tempfile1,
                                      execpath=args.execscythe,
                                      logfilepath=logfile1,
                                      quality=code_p1, min_keep=0,
                                      prior=args.scytheprior))
            mainlog.write("Scythe on fq1 complete. {!s}"
                          " reads with adapters removed.\n"
                          "Scythe fq1 completed at: {!s} seconds.\n".format(
                              stats.stats['scythe'][0], round(time() - time0, 2)))
######## ===== Scythe Adapter Trimming on Pair 2 ===================
        if paired_end and not args.skipscythe:
            tempfile2 = os.path.abspath("{!s}/scythetemp_p2_{!s}_{!s}"
                                        ".fastq".format(
                                            args.tempdir, timestamp, ndex))
            logfile2 = os.path.abspath("{!s}/scythe_p2_{!s}_{!s}.log".format(
                args.tempdir, timestamp, ndex))
            mainlog.write("Running Scythe on fq2...\n")
            stats.p2scythe(run_scythe(fq2, adapters=adapterpath, out=tempfile2,
                                      execpath=args.execscythe,
                                      logfilepath=logfile2,
                                      quality=code_p2, min_keep=0,
                                      prior=args.scytheprior))
            mainlog.write("Scythe on fq2 complete. {!s} reads with adapters"
                          "removed.\nScythe fq1 completed at"
                          " {!s} seconds.".format(stats.stats['scythe'][2],
                                                  round(time() - time0, 2)))
######## ===== Quality Trimming Setup ========================
        lowq_front = set([])
        lowq_end = set([])
        qual_offset = 32
        if code_p1 in ['solexa', 'illumina_1.3', 'illumina_1.5']:
            qual_offset = 64
        if trimqual[0]:
            lowq_front = set([chr(x) for x in range(qual_offset,
                                                    qual_offset + trimqual[0])])
        if trimqual[1]:
            lowq_end = set([chr(x) for x in range(qual_offset,
                                                  qual_offset + trimqual[1])])
######## ====== Open Files and Establish Trimming Counters =================
        if args.skipscythe:
            infq1 = open(fq1, 'r')
        else:
            infq1 = open(tempfile1, 'r')
        if paired_end:
            if args.skipscythe:
                infq2 = open(fq2, 'r')
            else:
                infq2 = open(tempfile2, 'r')
        while 1:
            read1_removed = ''
            read2_removed = ''
############ Read1 Filtering =============================================
            line1_header = infq1.readline().rstrip()
            line1_seq = infq1.readline().rstrip()
            infq1.readline()
            line1_qual = infq1.readline()
            if not line1_qual:
                break
            if args.cleanheader:
                line1_header = line1_header.split()[0]
            line1_qual = line1_qual.rstrip()
            len_read1 = len(line1_qual)
            stats.p1add('total', len_read1)
            if len_read1 < args.filterlength:
                read1_removed = 'tooshort'
            if not read1_removed and args.filterambig:
                if line1_seq.count('N') > args.filterambig:
                    read1_removed = 'ambig'
            if not read1_removed and args.filterquality:
                if sum([ord(x) - qual_offset for x in line1_qual
                        ])/(len_read1) < args.filterquality:
                    read1_removed = 'totalqual'
            if not read1_removed and args.filterlowinfo:
                mutinfo = mutual_info(line1_seq)
                if mutinfo >= args.filterlowinfo:
                    read1_removed = 'lowinfo'
############ Read2 Filtering ============================================
            if paired_end:
                line2_header = infq2.readline().rstrip()
                line2_seq = infq2.readline().rstrip()
                infq2.readline()
                line2_qual = infq2.readline()
                if not line2_qual:
                    break
                if args.cleanheader:
                    line2_header = line2_header.split()[0]
                line2_qual = line2_qual.rstrip()
                len_read2 = len(line2_qual)
                stats.p2add('total', len_read2)
                if read1_removed and args.filterunpaired:
                    read2_removed = 'unpaired'
                elif len_read2 < args.filterlength:
                    read2_removed = 'tooshort'
                if not read2_removed and args.filterambig:
                    if line2_seq.count('N') > args.filterambig:
                        read2_removed = 'ambig'
                if not read2_removed and args.filterquality:
                    if sum([ord(x) - qual_offset for x in line2_qual])/(
                            len_read2) < args.filterquality:
                        read2_removed = 'totalqual'
                if not read2_removed and args.filterlowinfo:
                    mutinfo = mutual_info(line2_seq)
                    if mutinfo >= args.filterlowinfo:
                        read2_removed = 'lowinfo'
            if read2_removed and (not read1_removed and args.filterunpaired):
                read1_removed = 'unpaired'
############### ============ Read1 trimming from 3' end ===============
            if not read1_removed:
                coord_end1 = len_read1 - trimfixed[1] - 1
                if args.trimpattern3:
                    for elem in args.trimpattern3:
                        if line1_seq.endswith(elem):
                            coord_end1 = len_read1 - (
                                max(len(elem), trimfixed[1]) - 1)
                            break
                for nuc, num in polynucs:
                    if len(line1_seq) <= num:
                        continue 
                    if line1_seq.endswith(nuc * num):
                        j = len(line1_seq) - num
                        while line1_seq[j] == nuc and j:
                            j -= 1
                        coord_end1 = j + 1
                while not read1_removed:
                    if coord_end1 <= args.filterlength + trimfixed[0]:
                        read1_removed = 'totalqual'
                        break
                    if line1_qual[coord_end1] not in lowq_end:
                        if not (set(line1_qual[coord_end1 - trimqualpad[1]:
                                               coord_end1:-1]) & lowq_end):
                            break
                    coord_end1 -= 1
################ Read1 trimming for 5' end ===============================
                coord_start1 = trimfixed[0] + 0
                for nuc, num in polynucs:
                    if len(line1_seq) <= num:
                        continue 
                    if line1_seq.startswith(nuc * num):
                        j = num + 0
                        while line1_seq[j] == nuc:
                            j += 1
                            if j >= len(line1_seq):
                                break
                        coord_start1 = j
                while not read1_removed:
                    if coord_start1 >= coord_end1:
                        read1_removed = 'totalqual'
                        break
                    if line1_qual[coord_start1] not in lowq_front:
                        if not (lowq_front & set(line1_qual[coord_start1:
                                                            coord_start1
                                                            + trimqualpad[0]])):
                            break
                    coord_start1 += 1
                if coord_end1 - coord_start1 < args.filterlength:
                    read1_removed = 'totalqual'
############ Read2 trimming of 3' end=============================
            if not read2_removed and paired_end:
                coord_end2 = len_read2 - trimfixed[1] - 1
                if args.trimpattern3:
                    for elem in args.trimpattern3:
                        if line2_seq.endswith(elem):
                            coord_end2 = len_read2 - (
                                max(len(elem), trimfixed[1]) - 1)
                            break
                for nuc, num in polynucs:
                    if len(line2_seq) <= num:
                        continue 
                    if line2_seq.endswith(nuc * num):
                        j = len(line2_seq) - num
                        while line2_seq[j] == nuc:
                            j -= 1
                            if j == 0:
                                break
                        coord_end2 = j + 1
                while not read2_removed:
                    if coord_end2 <= args.filterlength + trimfixed[0]:
                        read2_removed = 'totalqual'
                        break
                    if line2_qual[coord_end2] not in lowq_end:
                        if not (lowq_end &
                                set(line2_qual[coord_end2 - trimqualpad[1]:
                                               coord_end2])):
                            break
                    coord_end2 -= 1
################ Read2 trimming of 5' end ====================================
                coord_start2 = trimfixed[0] + 0
                for nuc, num in polynucs:
                    if len(line2_seq) <= num:
                        continue 
                    if line2_seq.startswith(nuc * num):
                        j = num + 0
                        while line2_seq[j] == nuc:
                            j += 1
                            if j >= len(line2_seq):
                                break
                        coord_start2 = j
                while not read2_removed:
                    if coord_start2 >= len_read2 or coord_start2 >= coord_end2:
                        read2_removed = 'totalqual'
                        break
                    if line2_qual[coord_start2] not in lowq_front:
                        if not (lowq_front & set(line2_qual[coord_start2:
                                                            coord_start2
                                                            + trimqualpad[0]])):
                            break
                    coord_start2 += 1

                if not read2_removed and (coord_end2 - coord_start2
                                          < args.filterlength):
                    read2_removed = 'tooshort'
            if read1_removed and (not read2_removed and args.filterunpaired):
                read2_removed = 'unpaired'
            if read2_removed and (not read1_removed and args.filterunpaired):
                read1_removed = 'unpaired'
            if not read1_removed:
                outfq1.write("{!s}\n{!s}\n+\n{!s}\n".format(
                    line1_header,
                    line1_seq[coord_start1:coord_end1 + 1],
                    line1_qual[coord_start1:coord_end1 + 1]))
                ntrim = len_read1 - (coord_end1 - coord_start1)
                if ntrim:
                    stats.p1add('trim', ntrim)
            else:
                stats.p1add(read1_removed, len_read1)
                if args.filterfile:
                    filterfq1.write("{!s}_{!s}\n{!s}\n+\n{!s}\n".format(
                        line1_header, read1_removed, line1_seq, line1_qual))
            if paired_end:
                if read2_removed:
                    stats.p2add(read2_removed, len_read2)
                    if args.filterfile:
                        filterfq2.write("{!s}_{!s}\n{!s}\n+\n{!s}\n".format(
                            line2_header, read2_removed, line2_seq, line2_qual))
                else:
                    outfq2.write("{!s}\n{!s}\n+\n{!s}\n".format(
                        line2_header,
                        line2_seq[coord_start2:coord_end2 + 1],
                        line2_qual[coord_start2:coord_end2 + 1]))
                    ntrim = len_read2 - (coord_end2 - coord_start2)
                    if ntrim:
                        stats.p2add('trim', ntrim)

        infq1.close()
        if paired_end:
            infq2.close()
        ndex += 1
######## FINISH ITERATION
    outfq1.close()
    if args.filterfile:
        filterfq1.close()
    if paired_end:
        outfq2.close()
        if args.filterfile:
            filterfq2.close()
#### ===== Write Statistics =======
    mainlog.write("\nfq1_stats\n{!s}\n\nfq2 stats\n{!s}\n".format(
        stats.stat_list(pair1=True), stats.stat_list(pair2=True)))
#### ==== Optional File Cleanup ======
    if args.cleanup != "none":
        os.remove(tempfile1)
        os.remove(tempfile2)
        if args.cleanup != "tempfastq":
            os.remove(logfile1)
            os.remove(logfile2)
            if args.cleanup != "exceptadapters":
                os.remove(adapterpath)
#### ===== Finish log and close ====
    mainlog.write("\nTotal running time: {!s} seconds.\n".format(
        round(time() - time0, 2)))
    mainlog.close()
    return ''

if __name__ == "__main__":
    main()
