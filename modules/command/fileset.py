import os
import pandas as pd
import sys
import shutil

class FileSet:
    def __init__(self):
        pass

    def run(self, main):
        self.create_transrate_paths(main)
        self.create_aligner_paths(main)
        self.set_file_names(main)

    def create_transrate_paths(self, main):
        # Create transrate directory
        if main.ASSEMBLY_RUN <= 1:
            for root, dirs, files in os.walk(main.OUTPUT):
                for file in files:
                    shutil.remove(os.path.join(root, file))
                for dir in dirs:
                    shutil.rmtree(os.path.join(root, dir))

        main.TRANSRATE_PATH = os.path.join(main.OUTPUT, 'transrate')
        main.ASSEMBLY_FILE = os.path.join(main.TRANSRATE_PATH, 'assembly.csv')
        main.LOG_PATH       = os.path.join(main.OUTPUT, 'logs')
        if not os.path.exists(main.TRANSRATE_PATH):
            os.makedirs(main.TRANSRATE_PATH)
        if not os.path.exists(main.LOG_PATH):
            os.makedirs(main.LOG_PATH)

    def create_aligner_paths(self, main):
        # Create aligner directory
        if not os.path.exists(main.OUTPUT):
            os.makedirs(main.OUTPUT)

        if main.ALIGNER == 'snap':
            main.ALIGNER_PATH  = os.path.join(main.OUTPUT, 'snap_' + main.ASSEMBLY_NAME)
            main.ALIGNER_INDEX = os.path.join(main.ALIGNER_PATH, 'index')
            main.BAM_ALIGNER   = os.path.join(main.ALIGNER_PATH, f'{main.ASSEMBLY_NAME}.bam')
            if not os.path.exists(main.ALIGNER_PATH):
                os.makedirs(main.ALIGNER_PATH)
            if not os.path.exists(main.ALIGNER_INDEX):
                os.makedirs(main.ALIGNER_INDEX)
        elif main.ALIGNER == 'star':
            main.ALIGNER_PATH  = os.path.join(main.OUTPUT, 'star_' + main.ASSEMBLY_NAME)
            main.ALIGNER_INDEX = os.path.join(main.ALIGNER_PATH, 'index')
            main.BAM_ALIGNER   = os.path.join(main.ALIGNER_PATH, f'{main.ASSEMBLY_NAME}_Aligned.out.bam')
            if not os.path.exists(main.ALIGNER_PATH):
                os.makedirs(main.ALIGNER_PATH)
            if not os.path.exists(main.ALIGNER_INDEX):
                os.makedirs(main.ALIGNER_INDEX)
        elif main.ALIGNER == 'bowtie2':
            main.ALIGNER_PATH  = os.path.join(main.OUTPUT, 'bowtie2_' + main.ASSEMBLY_NAME)
            main.ALIGNER_INDEX = os.path.join(main.ALIGNER_PATH, 'index')
            main.BAM_ALIGNER   = os.path.join(main.ALIGNER_PATH, f'{main.ASSEMBLY_NAME}.bam')
            if not os.path.exists(main.ALIGNER_PATH):
                os.makedirs(main.ALIGNER_PATH)
            if not os.path.exists(main.ALIGNER_INDEX):
                os.makedirs(main.ALIGNER_INDEX)
        else:
            # print('No aligner specified')
            # sys.exit(1)
            main.LOG.error_out(main, 'Aligner', 'No aligner specified')

        main.SALMON_PATH  = os.path.join(main.OUTPUT, 'salmon_' + main.ASSEMBLY_NAME)
        main.SALMON_QUANT = os.path.join(main.SALMON_PATH, 'quant.sf')
        main.BAM_SALMON   = os.path.join(main.SALMON_PATH, f'{main.ASSEMBLY_NAME}.postSample.bam')
        main.BAM_SORTED   = os.path.join(main.SALMON_PATH, f'{main.ASSEMBLY_NAME}.postSample.sorted.bam')
        main.TR_CSV       = os.path.join(main.TRANSRATE_PATH, f'{main.ASSEMBLY_NAME}.transrate.csv')
        if not os.path.exists(main.SALMON_PATH):
            os.makedirs(main.SALMON_PATH)

    def set_file_names(self, main):
        # Set contig file name based on assembly count
        # If multiple assemblies, use assembly name + '_contigs.csv'
        # If single assembly, use 'contigs.csv'
        main.CONTIG_FILE = os.path.join(main.TRANSRATE_PATH, 'contigs.csv') if not main.ASSEMBLY_MULTI else os.path.join(main.TRANSRATE_PATH, main.ASSEMBLY_NAME + '.contigs.csv')
        # if os.path.exists(main.CONTIG_FILE):
        #     os.remove(main.CONTIG_FILE)

        # Set assembly file name based on assembly count
        # if os.path.exists(self.ASSEMBLY_FILE) and not main.ASSEMBLY_MULTI:
        #     os.remove(self.ASSEMBLY_FILE)
        # If multiple assemblies, set assembly row to write into blank row of existing file
        if os.path.exists(main.ASSEMBLY_FILE) and main.ASSEMBLY_MULTI:
            main.ASSEMBLY_ROW = len(pd.read_csv(main.ASSEMBLY_FILE)) + 1
