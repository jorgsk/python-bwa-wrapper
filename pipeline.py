"""
Program for running bwa on nuclear reads from 5 different cell lines

"""
from IPython.Debugger import Tracer
debug = Tracer()

from subprocess import Popen, call

import os
import glob
# You tested having on read pair in the index and another not -- that works
# fine. You must also try the funny faste file output by bedtools.

def get_dsets(dset_dir, cell_lines):
    """
    Return the paths to the datsets for nuclear rna-seq
    """

    dsets = {}
    for line in glob.glob(dset_dir+'/*'):

    # I'm cheating for now!
    #for line in open('those_files', 'rb'):

        if 'NucleusPap' not in line or 'fastq' not in line:
            continue

        dname, fname = os.path.split(line.rstrip())

        for cl in cell_lines:

            if cl in fname:

                info = fname.split('Rd')[1]
                read_nr, more = info.split('Rep')
                biorep_nr = more[0]

                dset = cl+'_{0}'.format(biorep_nr)

                if dset in dsets:
                    dsets[dset].append(fname)
                else:
                    dsets[dset] = [fname]

    return dsets

def makedirs(dsets):
    """
    Make all output directories as needed
    """
    outp = 'output'
    if not os.path.isdir(outp):
        os.mkdir(outp)

    for dset in dsets:

        # dset
        dset_out = os.path.join(outp, dset)
        if not os.path.isdir(dset_out):
            os.mkdir(dset_out)

        # sai
        sai = os.path.join(dset_out, 'sai')
        if not os.path.isdir(sai):
            os.mkdir(sai)

        # sam
        sam = os.path.join(dset_out, 'sam')
        if not os.path.isdir(sam):
            os.mkdir(sam)


def bwa_launcher(fasta, cell_lines, dsets, dset_dir, my_bwa, speedrun):
    """
    # bwa index smallindex/smallref.fa

    # bwa requires two steps after indexing
    # step 1: make mapping tables

    # my options
    #  -I        The input is in the Illumina 1.3+ read format (quality equals ASCII-64).
    #small
    bwa aln -I smallindex/smallref.fa smallread1.fq > read1.sai
    bwa aln -I smallindex/smallref.fa smallread2.fq > read2.sai

    step 2: create bamfiles as output

    step 3: sort the  bamfile

    step 4: run rmdupliactes

    step 5: index the file

    You now have a file that is ready to be used for variant calling

    Put everything in a pipe to save disc-space and temporary files

    UPDATE now running rmdup without -S flag as it causes downstream errors
    """

    # prepare output directories
    makedirs(dsets)

    # make index if not made
    if not os.path.isfile(fasta+'.bwt'):
        # bwa will make the index in the correct dir
        cmd = [my_bwa, 'index', fasta]
        call(cmd) # generate index and wait until it's done

    # run bwa aln with 2 threads for each read-file

    for dset, fastqs in dsets.items():
        outp = os.path.join('output', dset)
        f0, f1 = fastqs

        f0 = os.path.join(dset_dir, f0)
        f1 = os.path.join(dset_dir, f1)

        # for testing!
        #f0 = '/users/rg/jskancke/phdproject/3UTR/UTR_exon_remapping/sample_data/k562Cyt1.fastq'
        #f1 = '/users/rg/jskancke/phdproject/3UTR/UTR_exon_remapping/sample_data/k562Cyt2.fastq'

        out0 = os.path.join(outp, 'sai', dset+'_0.sai')
        out1 = os.path.join(outp, 'sai', dset+'_1.sai')

        # skop the sai-generation if the files exist (hopefully they were
        # sucessfully generated)
        if not os.path.isfile(out0):
            out0handle = open(out0, 'wb')
            out1handle = open(out1, 'wb')

            cmd1 = [my_bwa, 'aln', '-t', '2', '-I', fasta, f0]
            cmd2 = [my_bwa, 'aln', '-t', '2', '-I', fasta, f1]

            print('Creating the sai files ...')
            p1 = Popen(cmd1, stdout=out0handle)
            p2 = Popen(cmd2, stdout=out1handle)

            p1.wait()
            p2.wait()

        out_bam = os.path.join(outp, 'sam', dset+'_nodup_sorted.bam')

        # NOTE pipe uncompressed -u sam for faster working
        # NOTE the bwa sampe includes all reads, even those who don't map.
        # You'll either have to ask it to output only those pairs that map or
        # you'll have to remove those pairs where both pairs don't map.
        # NOTE reemember -P for sampe, loads index into memory
        # NOTE for samtools sort, must specify -o for getting to stdout
        # To get unique hits, filter by...

        # If speedrun, grab only 1 mill reads
        if speedrun:

            cmd3 = '{6} sampe -P {0} {1} {2} {3} {4} | '\
                     'head -1000000 | '\
                     'samtools view -u -q 1 -F 0X4 -Sbh - | '\
                     'samtools sort -o - tmp_holder | '\
                     'samtools rmdup - {5}'\
                    .format(fasta, out0, out1, f0, f1, out_bam, my_bwa)
        else:

            cmd3 = '{6} sampe -P {0} {1} {2} {3} {4} | '\
                     'samtools view -u -q 1 -F 0X4 -Sbh - | '\
                     'samtools sort -o - tmp_holder | '\
                     'samtools rmdup - {5}'\
                    .format(fasta, out0, out1, f0, f1, out_bam, my_bwa)

        print('Filtering, sorting, and removing duplicates {0} ...'.format(dset))
        # run the sampe + samtools pipeline
        p = Popen(cmd3, shell=True)
        p.wait()

        # index the final product
        call(['samtools', 'index', out_bam])


def main():

    # where the fastq files are
    dset_dir = '/users/rg/projects/encode/scaling_up/whole_genome'\
            '/encode_DCC_mirror/wgEncodeCshlLongRnaSeq/releaseLatest'

    # the fasta-file of the 3UTR exons
    # NOTE: and how did you make this file? Using bed->fasta OK, but using what
    # bedfile?
    fasta = 'fasta/utr_exons.fasta'

    cell_lines = ['Gm12878', 'K562', 'Huvec', 'Helas3', 'Hepg2']
    #cell_lines = ['Gm12878', 'K562']
    #cell_lines = ['Huvec', 'Helas3', 'Hepg2']
    #cell_lines = ['Gm12878']

    dsets = get_dsets(dset_dir, cell_lines)

    #speedrun = True
    speedrun = False

    my_bwa = '/users/rg/jskancke/programs/other/bin/bwa'
    bwa_launcher(fasta, cell_lines, dsets, dset_dir, my_bwa, speedrun)


if __name__ == '__main__':
    main()
