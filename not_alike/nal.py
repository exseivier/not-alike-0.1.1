#!/usr/bin/env python3
import random
import click
#import utils as CMD
import not_alike.utils as CMD


@click.group()
def main():
    """
        Not-Alike: Command pipeline that identifies dissimilar regions of a target genome by comparing it to a genomes database.
    """
    pass

######################################
######          SEARCH          ######
######################################

@main.command()
@click.option('-g', '--genome', \
                help = 'Query genome FASTA file name', \
                required = True, \
                type = str)
@click.option('-ws', '--window-size', \
                help = 'Fragments from split genome will size equal to window size', \
                required = True, \
                type = int)
@click.option('-ss', '--step-size', \
                help = 'Step size in nucleotides the window takes between each cut along the genome', \
                required = True, \
                type = int)
@click.option('-db', '--database-file', \
                help = 'A name of a file that contains the pathway to BLAST_DB ver5 (*.db) files separated by return character', \
                required = True, \
                type = str)
@click.option('-e', '--evalue', \
                help = 'E-value cutoff', \
                required = True, \
                type = float)
@click.option('-i', '--identity', \
                help = 'Identity percentage cutoff', \
                required = True, \
                type = float)
@click.option('-q', '--qcov', \
                help = 'HSP query coverage cutoff', \
                required = True, \
                type = float)
@click.option('-t', '--task', \
                help = 'BLAST task [blastn | megablast | dc-megablast]', \
                required = True, \
                type = str)
@click.option('-nc', '--num-cores', \
                help = 'number of cores you want to use. It could be (PC total cores - 1)', \
                required = False, \
                default = 1,
                type = int)
@click.option('-c', '--comment', \
                help = 'Leave a comment enclosed by single quotes', \
                required = False, \
                default = 'Empty', \
                type = str)
@click.option('--quite-opposite', \
                is_flag = True, \
                help = 'Performs the opposite task of not-alike.', \
                type = bool)
def search(genome, window_size, step_size, database_file, evalue, identity, qcov, task, comment, quite_opposite, num_cores):
    """
        Searches for not alike fragments in query genome
    """
    PID = random.randrange(1, 9999999999)
    out_split = '.'.join(genome.split('.')[:-1]) + '_split_' + str(window_size) + '_' + str(step_size) + '.fasta'
    out_split = ''.join(out_split.split('/')[-1])
    input_split = 'input_split.' + str(PID) + '.fasta'
    print(out_split + ' was loaded.')
    print(input_split + ' was loaded.')

    CMD.check_path_exists('split_out')
    CMD.check_path_exists('blast_db')
    CMD.check_path_exists('blast_out')
    CMD.check_path_exists('ht2_idx')
    CMD.check_path_exists('mapping')
    CMD.check_path_exists('gtfs')
    CMD.check_path_exists('log')

    CMD.loggin_data(PID, genome, database_file, window_size, step_size, task, identity, qcov, evalue, comment)
    
    out_split = 'split_out/' + out_split
    input_split = 'split_out/' + input_split
    if not CMD.os.path.exists(out_split):
        CMD.split_genome(genome, window_size, step_size, out_split)
    else:
        print('A split-genome file was found!!!')

    CMD.copy_file(out_split, input_split)

    database_file_path = '/'.join(database_file.split('/')[:-1])

    db_files = CMD.load_lines(database_file)

    for f in db_files:
        dbf = '.'.join(f.split('.')[:-1]) + '.db'
        print('Blasting ' + dbf + ' ...')
        CMD.do_blast(input_split, database_file_path + '/' + dbf, 'blast_out/out.blast', evalue, identity, qcov, task, num_cores)
        print('Updating input_split')
        CMD.select_sequences(input_split, 'blast_out/out.blast', quite_opposite)

    print('BLASTn searching done!')
    print('Mapping on process.')
    
    CMD.mapping(genome, input_split)

    print('Assembly on process.')
    CMD.assembly(PID)

    print('Extracting sequences.')
    CMD.extseq(genome, PID)

    print('Doing assembly stats')
    CMD.do_assembly_stats('gtfs/nal_frags.' + str(PID) + '.fasta', PID)

    print('not-alike has finished.')

##############################################
######          DB MAKEBLAST            ######
##############################################

@main.command()
@click.option('-db', '--db-path', \
                help = 'Path to FASTA files database', \
                required = True, \
                type = str)
def db_makeblast(db_path):
    """
        Builds a BLAST_DB (version 5) database files
    """

    CMD.make_db(db_path)
    
##########################################
######          DB MAKEFILE         ######
##########################################

@main.command()
@click.option('-db', '--db-path', \
                help = 'Path to FASTA files database', \
                required = True, \
                type = str)
@click.option('-e', '--exclude', \
                help = 'A list of accession numbers from the organisms you want to exclude from database text file.', \
                required = False, \
                default = None, \
                type = str)
@click.option('-i', '--include', \
                help = 'A list of accession numbers from the organism you want to include in database text file.', \
                required = False, \
                default = None, \
                type = str)
@click.option('-o', '--out', \
                help = 'Output file name', \
                required = True, \
                type = str)
def db_makefile(db_path, exclude, include, out):
    """
        Creates the database text file which contains the BLAST_DB files paths.
    """
    if exclude == None and include == None:
        print('ERROR, Exclude and include options are empty')
        return 1
    if exclude != None and include != None:
        print('ERROR, Only one option between exclude or include is allowed')
        return 1
    if exclude != None:
        exclude = exclude.split(',')
    if include != None:
        include = include.split(',')
    CMD.make_txtfiledb(db_path, exclude, include, out)

######################################
######          SHOW DB         ######
######################################

@main.command()
@click.option('-p', '--db-path', \
                help = 'Genomes database path', \
                required = True, \
                type = str)
def show_db(db_path):
    """
        Shows metadata of database [accession number, organism name and organism taxon id]
    """

    assembly_report = db_path + '/assembly_data_report.jsonl'
    assembly_report_tsv = db_path + '/assembly_data_report.tsv'
    if not CMD.file_exists(assembly_report_tsv):
        CMD.dataformat_tsv(assembly_report, assembly_report_tsv)

    print(assembly_report_tsv)

    CMD.print_table(assembly_report_tsv)

##########################################
######          SHOW EXP            ######
##########################################

@main.command()
@click.option('--sort-by', \
                help = 'Criteria to sort values.', \
                required = False, \
                default = 'DATE', \
                type = str)
def show_exp(sort_by):
    """
        Shows information about epxeriments stored in the current working directory
    """
    CMD.show_exp_info(sort_by)

##########################################
######          ASSM STATS          ######
##########################################

@main.command()
@click.option('-f', '--file-name', \
                help = 'FASTA file name', \
                required = True, \
                type = str)
@click.option('-pid', \
                help = 'Process ID', \
                required = False, \
                default = 00, \
                type = int)
def assm_stats(file_name, pid):
    """
        Calculates assembly statistics such as: Mean, Median, N50 and L50.
    """
    CMD.do_assembly_stats(file_name, pid)

##############################################
######          PRIMER SELECT           ######
##############################################

@main.command()
@click.option('--input-file', \
                help = 'Input fasta file', \
                required = True, \
                type = str)
@click.option('--opt-size', \
                help = 'Optimum primer size (nt)', \
                required = True, \
                type = int)
@click.option('--opt-gc', \
                help = 'Optimum GC percentage (%)', \
                required = True, \
                type = float)
@click.option('--opt-tm', \
                help = 'Optimum melting point (°C)', \
                required = True, \
                type = float)
@click.option('--product-size', \
                help = 'Expected product zise (bp) [i.e. 75-100]', \
                required = True, \
                type = str)
@click.option('--template-size-range', \
                help = 'Template sequence size range (bp) [i.e. 750-1000]', \
                required = True, \
                type = str)
def search_primers(input_file, opt_size, opt_gc, opt_tm, product_size, template_size_range):
    '''
        Selects the best fitted primer sequences based on input arguments.
    '''
    CMD.find_primers(input_file, opt_size, opt_gc, opt_tm, product_size, template_size_range)

if __name__ == '__main__':
    main()
