import os
import subprocess as sup
import pandas as pd


########################################
######                            ######
######      HELPER FUNCTIONS      ######
######                            ######
########################################


######      HELPER FUNCTION FOR SPLIT-GENOME

def assert_directory(path):
    """
        Asserts directory exists
    """
    path = path.split('/')
    file_exists = os.path.exists('/'.join(path))
    if file_exists:
        print('File exists')
    else:
        print('File does not exists!!!!')
        print('Creating directory and / or file')
        if len(path) > 1:
            directory = ''.join(path[:-1])
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                pass
            path = '/'.join([directory, path[-1]])
            with open(path, 'w') as fh:
                fh.close()
        else:
            with open(''.join(path), 'w') as fh:
                fh.close()

def __load_seqs(filename):
    """

    """
    FHIN = open(filename, "r")

    line = FHIN.readline()
    seqs = {}
    while line:
        line = line.strip("\n")
        if line[0] == ">":
            header = line
            seqs[header] = []
        else:
            seqs[header].append(line)

        line = FHIN.readline()
    out_seqs = {}
    for head, seq in seqs.items():
        out_seqs[head] = "".join(seq)

    FHIN.close()
    seqs = None
    return out_seqs


def load_lines(infile):
    """
        Loads lines from file and stores them in a list
    """
    lines = []
    with open(infile, 'r') as fh:
        for line in fh:
            lines.append(line.strip())
        
        fh.close()

    return lines

def __split(seqs, size, step):
    """

    """
    out_seqs = {}
    seq_counter = 0
    for head, seq in seqs.items():
        lseq = len(seq)
        i = 0
        start = i
        while start < lseq:
            end = start + size - 1
            out_seqs[">fragment_" + str(seq_counter)] = seq[start:end]
            start = start + step
            i += 1
            seq_counter += 1

    return out_seqs

def __write_seqs(seqs, outfile):
    """

    """
    FHOUT = open(outfile, "w+")
    
    big_string = ""
    for head, seq in seqs.items():
        big_string += head + "\n" + seq + "\n"

    FHOUT.write(big_string)
    FHOUT.close()

##################################################

######      HELPER FUNCTION FOR SELECT-SEQUENCES

def __load_headers(hd_file):
    """
        Load the headers of the sequences that hitted a subject in genomes database.
    """
    
    FHIN = open(hd_file, "r")

    line = FHIN.readline()
    heads = []
    while line:

        line = line.strip("\n")
        heads.append(line)
        line = FHIN.readline()

    FHIN.close()
    return heads


def __select_seqs(seqs, heads, quite_opposite):
    """
        Selects those sequences which header is not in heads list.
    """

    lsq = len(seqs)
    print(str(lsq) + ' current sequences')
    heads = list(set(heads))
    lhd = len(heads)
    if quite_opposite:
        print(str(lhd) + ' similar sequences')
    else:
        print(str(lhd) + ' dropped secuences')
    if lhd <= 0:
        return seqs
    list_seqs_keys = [x[1:] for x in list(seqs.keys())]
    list_seqs_keys = list(set(list_seqs_keys))

    heads_seqs_keys = list_seqs_keys + heads
    heads_seqs_keys = sorted(heads_seqs_keys)
    
    pivot = 0
    nextp = pivot + 1

    # Modified section for version 0.1.1

    while pivot < len(heads_seqs_keys)-1:
        if heads_seqs_keys[pivot] == heads_seqs_keys[nextp]:
            if quite_opposite:
                heads_seqs_keys.pop(pivot)
                pivot += 1
                nextp = pivot + 1
            else:
                heads_seqs_keys.pop(pivot)
                heads_seqs_keys.pop(pivot)
        else:
            if quite_opposite:
                heads_seqs_keys.pop(pivot)
            else:
                pivot += 1
                nextp = pivot + 1
    

    out_seqs = {}

    for head in heads_seqs_keys:
        head = '>' + head
        out_seqs[head] = seqs[head]

    #####################################

    print(str(len(out_seqs)) + ' maintained sequences')
    return out_seqs

#####################################################################

def file_exists(infile):
    """
        Checks if path is a regular file
    """

    if os.path.exists(infile) and os.path.isfile(infile):
        return True
    else:
        return False

def check_path_exists(path):
    """
        Checks path exist
    """
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        print(path + ' folder exists!')


def copy_file(source, destiny):
    """
        Copies a file
    """
    p = sup.Popen(['cp', source, destiny], stdout = sup.PIPE, stderr = sup.PIPE)
    p.communicate()
    p.kill()
    if os.path.exists(destiny):
        print(source + ' was copied to ' + destiny)
    else:
        print('Copy process | Error!!!')

def rm_file(path):
    """
        Removes the file specifyied in path
    """
    if os.path.exists(path):
        os.remove(path)
    else:
        print(path + ' not found.')


######      HELPER FUNCTIONS FOR MAPPPING


def __ht2idx_ready():
    """
        Checks if Hisat2 database is allready formated.
    """

    for fname in os.listdir('ht2_idx'):
        if fname.endswith('.1.ht2'):
            return True
        
    return False


def __index(ref_genome, fname_suffix):
    """
        Prepare an index of reference genome (a.k.a. query genome)
    """

    p = sup.Popen(['hisat2-build', \
                ref_genome, \
                'ht2_idx/' + fname_suffix], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()


def __map(ref_genome, input_split, fname_suffix):
    """
        Maps sequences to reference genome
    """

    pid = input_split.split('.')[-2]
    mapping_out_sam = 'mapping/nal_frags.' + pid + '.sam'
    mapping_out_bam = 'mapping/nal_frags.' + pid + '.bam'
    mapping_out_sbam = 'mapping/nal_frags.' + pid + '.sort.bam'
    p = sup.Popen(['hisat2', \
                '-x', 'ht2_idx/' + fname_suffix, \
                '--no-temp-splicesite', \
                '--no-spliced-alignment', \
                '-f', \
                '-U', input_split, \
                '-S', mapping_out_sam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)

    p.communicate()
    p.kill()
    p = sup.Popen(['samtools', 'view', \
                '-b', '-h', \
                '-f', '0x0', '-F', '0x100', \
                '-o', mapping_out_bam, \
                mapping_out_sam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()
    p = sup.Popen(['samtools', 'sort', \
                '-o', mapping_out_sbam, \
                '-O', 'BAM', \
                '--reference', ref_genome, \
                mapping_out_bam], \
                stdout = sup.PIPE, \
                stderr = sup.PIPE)
    p.communicate()
    p.kill()

#########################################################

######      HELPER FUNCTION FOR ASSEMBLY

def __do_assembly(pid):
    p = sup.Popen(['stringtie', 'mapping/nal_frags.' + str(pid) + '.sort.bam', \
                    '-L', '-s', '1', '-g', '50', \
                    '-o', 'gtfs/nal_frags.' + str(pid) + '.gtf'], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)

    p.communicate()
    p.kill()


##########################################################


######      HELPER FUNCTION FOR EXTSEQ

def __extract_sequences(genome, pid):
    tmp_gff_file = 'gtfs/tmp.gff'
    fasta_out = 'gtfs/nal_frags.' + str(pid) + '.fasta'
    p = sup.Popen(['gffread', 'gtfs/nal_frags.' + str(pid) + '.gtf', \
                    '-o', tmp_gff_file], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)
    p.communicate()
    p.kill()

    p = sup.Popen(['gffread', tmp_gff_file, \
                    '-g', genome, \
                    '-w', fasta_out], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)
    p.communicate()
    p.kill()



######################################
######                          ######
######      MAIN FUNCTIONS      ######
######                          ######
######################################


def split_genome(in_file, size, step_size, out_file):
    """
        Split query genome in fragments of determined size and at each determined step
    """
    print("Loading genome...")
    seqs = __load_seqs(in_file)
    print("Splitting genome...")
    seqs = __split(seqs, size, step_size)
    print("Writting to file...")
    __write_seqs(seqs, out_file)

    print(f"Spliting genome {in_file}")


def do_blast(query, db_file, out_blast, evalue, idt, qcov, task):
    """
        Performs a BLASTn task.
    """
    
    p = sup.Popen(['blastn', \
                    '-query', query, \
                    '-db', db_file, \
                    '-out', out_blast, \
                    '-outfmt', '6 qseqid', \
                    '-task', task, \
                    '-perc_identity', str(idt), \
                    '-qcov_hsp_perc', str(qcov), \
                    '-evalue', str(evalue), \
                    '-max_target_seqs', str(1)], \
                    stdout = sup.PIPE, \
                    stderr = sup.PIPE)

    p.communicate()
    p.kill()

def select_sequences(in_file, hd_file, quite_opposite):
    """
        Selects those sequences that hit a subject in blast searching
    """
    
    seqs = __load_seqs(in_file)
    
    heads = __load_headers(hd_file)

    seqs = __select_seqs(seqs, heads, quite_opposite)

    __write_seqs(seqs, in_file)


def mapping(ref_genome, input_split):
    """

    """
    fname_suffix = ref_genome.split('/')[-1]
    fname_suffix = '.'.join(fname_suffix.split('.')[:-1])

    if not __ht2idx_ready():
        print('ht2 index not found.')
        print('Indexing.')
        __index(ref_genome, fname_suffix)
    else:
        print('ht2 index found.')
    
    print('Mapping to reference genome.')
    __map(ref_genome, input_split, fname_suffix)
    print('Mapping finished.')



def assembly(pid):
    """
        Assembles fragments of query genome using the genome-guided procedure.
    """

    __do_assembly(pid)

def extseq(genome, PID):
    """
        Extracts sequences from genome using gff or gtf input file.
    """
    __extract_sequences(genome, PID)

def dataformat_tsv(assembly_report, assembly_report_tsv):
    """
        Transforms jsonl to tsv format from ncbi dataset
    """
#    print(assembly_report)
#    print(assembly_report_tsv)
    fields = 'accession,organism-name,organism-tax-id'
    with open(assembly_report_tsv, 'w') as FHOUT:
        p = sup.Popen(['dataformat', 'tsv', 'genome', \
                        '--inputfile', assembly_report, \
                        '--fields', fields, \
                        '--force'],
                        stdout = FHOUT,
                        stderr = sup.PIPE)
        out, err = p.communicate()
#        print('hello world! ' + str(len(err)))
        p.kill()
        FHOUT.close()

def make_db(db_path):
    """
        Formats FASTA files to BLAST_DB format
    """
    
    dataset_catalog = db_path + '/dataset_catalog.json'
    assemblies = pd.read_json(dataset_catalog)
    assemblies = assemblies['assemblies']
    for assembly in assemblies:
        for fastafile in assembly['files']:
            if fastafile['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
                inputfile = db_path + '/' + fastafile['filePath']
                outputfile = '.'.join(inputfile.split('.')[:-1]) + '.db'
                p = sup.Popen(['makeblastdb', \
                            '-in', inputfile, \
                            '-dbtype', 'nucl', \
                            '-parse_seqids', \
                            '-out', outputfile],
                            stdout = sup.PIPE,
                            stderr = sup.PIPE)
                stdout, stderr = p.communicate()
                print(stdout)
                print(stderr)
                p.kill()

def make_txtfiledb(db_path, exclude, out):
    """
        Creates the text file which contains the path to BLAST_DB files
        you want them to be in the database.
        exclude is a list with the accession number of the organisms
        you want them to be excluded from database.
        db_path is the path where BLAST_DB files are.
    """
    
    catalog_file = db_path + '/dataset_catalog.json'
    output_file = db_path + '/' + out
    catalog = pd.read_json(catalog_file)
    assemblies = catalog['assemblies']
    out_files = []
    for assembly in assemblies:
        if 'accession' in assembly.keys():
            if assembly['accession'] in exclude:
                continue
            else:
                for _file in assembly['files']:
                    if _file['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
                        out_files.append(_file['filePath'])
                    else:
                        continue
        else:
            continue

    with open(output_file, 'a') as FHIN:
        FHIN.write('\n'.join(out_files))
        FHIN.close()





def loggin_data(pid, genome, database_file, window_size, step_size, task, identity, qcov, evalue, comment):
    """
        Stores experiment information in a log.
    """
    date = pd.to_datetime('today')
    date = '-'.join([str(date.year), str(date.month), str(date.day)])
    if not file_exists('log/experiments.log'):
        with open('log/experiments.log', 'a') as FHIN:
            FHIN.write('PID\tQSEQ\tDB_FILE\tWS\tSS\tTASK\tIDENT\tQCOV\tEVAL\tCOMMENT\tDATE\n')
            FHIN.write(str(pid) + '\t' \
                        + genome + '\t' \
                        + database_file + '\t' \
                        + str(window_size) + '\t' \
                        + str(step_size) + '\t' \
                        + task + '\t' \
                        + str(identity) + '\t' \
                        + str(qcov) + '\t' \
                        + str(evalue) + '\t' \
                        + comment + '\t' \
                        + date + '\n')
            FHIN.close()
    else:
        with open('log/experiments.log', 'a') as FHIN:
            FHIN.write(str(pid) + '\t' \
                        + genome + '\t' \
                        + database_file + '\t' \
                        + str(window_size) + '\t' \
                        + str(step_size) + '\t' \
                        + task + '\t' \
                        + str(identity) + '\t' \
                        + str(qcov) + '\t' \
                        + str(evalue) + '\t' \
                        + comment + '\t' \
                        + date + '\n')
            FHIN.close()

def print_table(tsv_file):
    """
        Prints on screen a pandas data frame.
    """

    tbl = pd.read_table(tsv_file, sep = '\t', header = None)
    print(tbl)

def show_exp_info(sort_by):
    """
        Prints on screen log/experiments.log file content.
    """
    sort_by = sort_by.split(',')   
    if os.path.exists('log/'):
        tbl = pd.read_table('log/experiments.log', sep = '\t')
        tbl.isetitem(10, [pd.Timestamp(x) for x in tbl.loc[:, 'DATE']])
        print(tbl.sort_values(by = sort_by))
    else:
        print('Unable to find log/ folder.')

