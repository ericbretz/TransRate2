import os
import sys
from pprint                     import pprint
from pathlib                    import Path
from core.utils.printout        import PrintOut
from core.alignments.bowtie2    import Bowtie2
from core.alignments.hisat2     import Hisat2
from core.alignments.salmon     import Salmon
from core.alignments.samtools   import Samtools
from core.contighub             import ContigHub
from core.assemblyhub           import AssemblyHub
from core.csvProcess            import CSV
from core.utils.logging        import TransRateLogger
import warnings
warnings.filterwarnings("ignore")
import pandas as pd

class DataHub:
    def __init__(self, args):
        # TransRate2
        self.threads            = args.threads                                                       # Number of threads
        self.mode_multi         = args.mode_multi                                                    # Multi-mode flag
        self.mode               = args.mode                                                          # 0: assembly only
                                                                                                     # 1: assembly and single-end reads
                                                                                                     # 2: assembly and paired-end reads
                                                                                                     # 3: BAM + assembly
                                                                                                     # 4: BAM only
        # Input files
        self.assembly_file      = args.assembly                                                      # Assembly file(s)
        self.assembly_name      = os.path.basename(args.assembly) if args.assembly else None         # Assembly name for CSV
        self.assembly_base_name = os.path.splitext(self.assembly_name)[0] if self.assembly_name \
                                  else os.path.splitext(os.path.basename(args.bam))[0] if args.bam \
                                  else 'analysis'                                                    # Assembly base name without extension
        self.left_reads         = args.left                                                          # Left reads file
        self.right_reads        = args.right                                                         # Right reads file
        self.reference_file     = args.reference                                                     # Reference file
        self.bam_file           = args.bam                                                           # BAM file

        # Aligner settings
        self.aligner_name       = args.aligner                                                       # Alignment tool (bowtie2/hisat2)
        self.bowtie2            = args.bowtie2                                                       # Bowtie2 flag
        self.hisat2             = args.hisat2                                                        # Hisat2 flag

        # Terminal settings
        self.log                = getattr(args, 'log', 3)                                            # Log level
        self.quiet              = getattr(args, 'quiet', False)                                      # Quiet mode
        self.nocolor            = getattr(args, 'nocolor', False)                                    # No color mode

        # Developer flags
        self.clutter            = args.clutter                                                       # Clutter flag
        self.debug              = getattr(args, 'debug', False)                                      # Debug flag

        # Initialize printout with colors and quiet mode
        highlight_color         = getattr(args, 'highlight_color', '\033[33m')
        background_color        = getattr(args, 'background_color', '\033[43m')
        self.printClass         = PrintOut(args.log, highlight_color, background_color)              # Printout class
        self.printClass.set_quiet(self.quiet)                                                        # Set quiet mode
        self.printClass.set_nocolor(self.nocolor)                                                    # Set nocolor mode
        self.printout           = self.printClass.printout                                           # Printout function

        # Master dictionary
        self.dict_all                  = {}                                                          # Master dictionary
        self.dict_dir                  = self._setup_directories(args)                               # Directory dictionary
        self.dict_file                 = self._setup_files()                                         # File dictionary
        self.dict_info                 = self._setup_info()                                          # Info dictionary
        self.contigDF  :pd.DataFrame   = pd.DataFrame()
        self.assemblyDF:pd.DataFrame   = pd.DataFrame()

        # Headers
        self.cHeaders = self._setup_cHeaders()                                                       # Contig headers
        self.aHeaders = self._setup_aHeaders()                                                       # Assembly headers

        # DataFrames
        self._initialize_dataframes()

        # Contig stuff
        self.basesDct  = {}
        self.refCount  = 0
        self.refList   = []
        self.readCount = 0

        # Multi-assembly
        self.multi_assembly = False

        # Logging system
        self.transrate_logger = TransRateLogger(
            log_dir       = self.dict_dir['logs'],
            assembly_name = self.assembly_base_name,
            quiet         = self.quiet
        )

        # Log setup
        self.transrate_logger.log_stage_start("TransRate2 Analysis", {
            "assembly": self.assembly_name,
            "mode"    : self.mode,
            "aligner" : self.aligner_name,
            "threads" : self.threads
        })

        # pprint(self.__dict__)

    def _setup_directories(self, args):
        assembly_name  = self.assembly_name.split('.')[0] if self.assembly_name else 'assembly'
        base_output    = Path(args.output_dir) if args.output_dir else Path.cwd()
        transrate2_dir = base_output.joinpath('TransRate2')

        base = {
            'input'             : Path.cwd(),
            'output'            : base_output,
            'transrate2'        : transrate2_dir,
            'results'           : transrate2_dir.joinpath('results').joinpath(assembly_name) if self.mode_multi else transrate2_dir.joinpath('results'),
            'logs'              : transrate2_dir.joinpath('logs').joinpath(assembly_name) if self.mode_multi else transrate2_dir.joinpath('logs'),
            'temp'              : transrate2_dir.joinpath('temp').joinpath(assembly_name) if self.mode_multi else transrate2_dir.joinpath('temp'),
            }

        temp_dirs = {
            'temp_aligner'      : base['temp'].joinpath('alignments').joinpath(self.aligner_name),
            'temp_aligner_index': base['temp'].joinpath('alignments').joinpath(self.aligner_name).joinpath('index'),
            'temp_salmon'       : base['temp'].joinpath('salmon'),
            'temp_bam'          : base['temp'].joinpath('bam'),
            'temp_analysis'     : base['temp'].joinpath('analysis'),
            'temp_reference'    : base['temp'].joinpath('reference'),
        } if self.mode > 0 else {}

        final = {**base, **temp_dirs}
        for key, path in final.items():
            path.mkdir(parents=True, exist_ok=True)
        return final

    def _setup_files(self):
        files = {
            'assembly'          : self.assembly_file,
            'assembly_name'     : self.assembly_name,
            'assembly_base'     : self.assembly_base_name,
            'left'              : self.left_reads,
            'right'             : self.right_reads,
            'single'            : self.left_reads if self.left_reads else self.right_reads,
            'reference'         : self.reference_file,
            'bam'               : self.bam_file,
        }
        aligner = {}
        if self.mode > 0:
            if self.bam_file:
                aligner = {
                    'aligner_bam'       : self.bam_file,
                    'salmon_bam'        : os.path.join(self.dict_dir['temp_salmon'], self.assembly_base_name + '_' + 'postSample.bam'),
                    'samtools_bam'      : os.path.join(self.dict_dir['temp_bam'], self.assembly_base_name + '_sorted.bam'),
                }
            else:
                aligner = {
                    'aligner_prefix'    : os.path.join(self.dict_dir['temp_aligner_index'], self.assembly_base_name),
                    'aligner_bam'       : os.path.join(self.dict_dir['temp_aligner'], self.assembly_base_name + '_aligned.bam') \
                                          if self.aligner_name == 'bowtie2' else os.path.join(self.dict_dir['temp_aligner'], self.assembly_base_name + '_aligned.sam'),
                    'salmon_bam'        : os.path.join(self.dict_dir['temp_salmon'], self.assembly_base_name + '_' + 'postSample.bam'),
                    'samtools_bam'      : os.path.join(self.dict_dir['temp_bam'], self.assembly_base_name + '_sorted.bam'),
                }

        final = {**files, **aligner}
        return final

    def _setup_info(self):
        info = {
            'mode'          : self.mode,
            'threads'       : self.threads,
            'aligner'       : self.aligner_name,
            'multi_mode'    : self.mode_multi,
            'clutter'       : self.clutter,
            'quiet'         : self.quiet,
            'debug'         : self.debug,
        }
        return info

    def _setup_cHeaders(self):
        base = [
            'name'        , 'length'     , 'gcCount'     , 'pGC'            , 'orfLength'    ,
            'fragments'   , 'softclipped', 'pSoftclipped', 'basesUncovered' , 'pBasesCovered', 'pSeqTrue', 
            'effLength'   , 'effCount'   , 'tpm'         , 'coverage'       , 'sCnuc'        , 'sCcov'
        ]

        paired = [
            'bridges'      , 'properPair', 'good'      , 'pGood'     , 'score',
            'pNotSegmented', 'sCord'     , 'sCseg'     , 'bothMapped'
        ]

        combined_headers = base + paired if self.mode == 2 else base

        final_order = [
            'name'         , 'length'       , 'fragments' , 'gcCount'    , 'pGC'         , 'basesUncovered',
            'pBasesCovered', 'bridges'      , 'bothMapped', 'properPair' , 'good'        , 'pGood'         ,
            'orfLength'    , 'pNotSegmented', 'pSeqTrue'  , 'softclipped', 'pSoftclipped', 'effLength'     ,
            'effCount'     , 'tpm'          , 'coverage'  , 'sCnuc'      , 'sCcov'       , 'sCord'         ,
            'sCseg'        , 'score'
        ]
        cHeaders = [header for header in final_order if header in combined_headers]
        return cHeaders

    def _setup_aHeaders(self):
        base = [
                'assembly'  , 'nSeqs'         , 'bases'    , 'smallest' , 'largest', 'basesN'  ,
                'meanLength', 'medianLength'  , 'stdLength', 'nUnder200', 'nOver1k', 'nOver10k',
                'nWithOrf'  , 'meanOrfPercent', 'n90'      , 'n70'      , 'n50'    , 'n30'     ,
                'n10'       , 'gcCount'       , 'pGC'      , 'basesN'   , 'pN'     
                ]
        single = [
                'fragments'         , 'fragmentsMapped'  , 'pFragmentsMapped' , 'softclipped'      , 
                'pSoftclipped'      , 'basesUncovered'   , 'pBasesUncovered'  , 'contigsUncovBase' ,
                'pContigsUncovbase' , 'contigsUncovered' , 'pContigsUncovered', 'contigsLowcovered', 
                'pContigsLowcovered',
                
        ]

        paired = [
                'goodMappings'    , 'bothMapped'       , 'pGoodMappings', 'badMappings' , 'potentialBridges', 
                'contigsSegmented', 'pContigsSegmented', 'score'        , 'optimalScore', 'cutoff'          , 
                'weighted'        , 'goodContigs'      , 'badContigs'   ,
        ]

        reference = [
                'CRBBhits'       , 'nContigsWithCRBB', 'pContigsWithCRBB', 'nRefsWithCRBB'    , 'pRefsWithCRBB',
                'rbhPerReference', 'nRefsWithCRBB'   , 'pRefsWithCRBB'   , 'cov25'            , 'pCov25'       ,
                'cov50'          , 'pCov50'          , 'cov75'           , 'pCov75'           , 'cov85'        , 
                'pCov85'         , 'cov95'           , 'pCov95'          , 'referenceCoverage'
        ]

        final_order = [
                'assembly'         , 'nSeqs'            , 'bases'            , 'smallest'          , 'largest',
                'meanLength'       , 'medianLength'     , 'stdLength'        , 'nUnder200'         , 'nOver1k'           ,
                'nOver10k'         , 'nWithOrf'         , 'meanOrfPercent'   , 'n90'               , 'n70'               ,
                'n50'              , 'n30'              , 'n10'              , 'gcCount'           , 'pGC'               ,
                'basesN'           , 'pN'               , 'fragments'        , 'fragmentsMapped'   , 'bothMapped'        ,
                'pFragmentsMapped' , 'softclipped'      , 'pSoftclipped'     , 'goodMappings'      , 'pGoodMappings'     ,
                'badMappings'      , 'potentialBridges' , 'basesUncovered'   , 'pBasesUncovered'   , 'contigsUncovBase'  ,
                'pContigsUncovbase', 'contigsUncovered' , 'pContigsUncovered', 'contigsLowcovered' , 'pContigsLowcovered', 
                'contigsSegmented' , 'pContigsSegmented', 'goodContigs'      , 'badContigs'        , 'cutoff'            ,
                'weighted'         , 'score'            , 'optimalScore'     , 'CRBBhits'          , 'nContigsWithCRBB'  ,
                'pContigsWithCRBB' , 'nRefsWithCRBB'    , 'pRefsWithCRBB'    , 'cov25'             , 'pCov25'            , 
                'cov50'            , 'pCov50'           , 'cov75'            , 'pCov75'            , 'cov85'             ,
                'pCov85'           , 'cov95'            , 'pCov95'           , 'referenceCoverage'
        ]

        if self.mode == 2:
            combined_headers = base + single + paired
        elif self.mode == 1:
            combined_headers = base + single
        else:
            combined_headers = base

        if self.reference_file:
            combined_headers = combined_headers + reference

        aHeaders = [header for header in final_order if header in combined_headers]
        return aHeaders

    def _initialize_dataframes(self):
        if self.mode > 0:  
            self.setup_contig_dataframe()

        self.assemblyDF = pd.DataFrame(columns=self.aHeaders)

        if self.assembly_name:
            assembly_val = {'assembly': self.assembly_base_name}
            self.assemblyDF.loc[0] = assembly_val

    def get_cHeaders(self):
        return self.cHeaders

    def get_aHeaders(self):
        return self.aHeaders

    def setup_contig_dataframe(self):
        if self.mode > 0 and self.assembly_file:
            try:
                from pysam import FastaFile
                fasta = FastaFile(self.assembly_file)
                refs = fasta.references

                contig_data = []
                for ref in refs:
                    contig_row = {'name': ref, 'length': fasta.get_reference_length(ref)}
                    for header in self.cHeaders:
                        if header not in contig_row:
                            contig_row[header] = None
                    contig_data.append(contig_row)

                self.contigDF = pd.DataFrame(contig_data, columns=self.cHeaders)
                fasta.close()

            except ImportError:
                self.printout('warning', 'pysam not available for contig DataFrame setup')
            except Exception as e:
                self.printout('warning', f'Error setting up contig DataFrame: {e}')

    def run(self):
        try:
            if self.mode == 0:          # Assembly only
                self.transrate_logger.log_stage_start("Assembly Analysis")
                self.assembly_run()
            elif self.mode == 4:        # BAM-only mode (no assembly)
                self.transrate_logger.log_stage_start("BAM-Only Processing")
                self.bam_run()
                self.bam_analysis_run()
            elif self.mode > 0:         # Assembly with reads (modes 1, 2, 3)
                if self.bam_file:       # BAM provided
                    self.transrate_logger.log_stage_start("BAM Processing")
                    self.bam_run()
                else:
                    self.transrate_logger.log_stage_start("Read Alignment", {"aligner": self.aligner_name})
                    self.aligner_run()

                self.salmon_run()
                self.samtools_run()
                self.contig_run()
                self.assembly_run()

            if self.reference_file:     # Reference
                self.reference_run()

            self.transrate_logger.log_stage_start("CSV Processing")
            self.process_csv_files()

            if self.clutter:
                self._cleanup_temp_files()

        finally:
            self.transrate_logger.finalize()

    def aligner_run(self):
        if self.aligner_name == 'bowtie2':
            self.bowtie2 = Bowtie2(self.dict_dir, self.dict_file, self.dict_info, self.printClass, self.transrate_logger)
            self.bowtie2.index()
            self.bowtie2.align()
        elif self.aligner_name == 'hisat2':
            self.hisat2 = Hisat2(self.dict_dir, self.dict_file, self.dict_info, self.printClass, self.transrate_logger)
            self.hisat2.index()
            self.hisat2.align()

    def salmon_run(self):
        self.salmon = Salmon(self.dict_dir, self.dict_file, self.dict_info, self.printClass, self.transrate_logger)
        self.salmon.quant()

    def samtools_run(self):
        self.samtools = Samtools(self.dict_dir, self.dict_file, self.dict_info, self.printClass, self.transrate_logger)
        self.samtools.sort()
        self.samtools.index()

    def bam_run(self):
        self.printout('subtitle', 'BAM Processing')
        self._validate_bam_file()

    def _validate_bam_file(self):
        if not os.path.exists(self.bam_file):
            self.printout('error', f'BAM file does not exist: {self.bam_file}')
            sys.exit(1)

        bam_size = self._get_bam_file_size()
        self.transrate_logger.log_progress(f"Using provided BAM file: {self.bam_file} ({bam_size})")
        self.printout('info', f'Using BAM file: {os.path.basename(self.bam_file)} ({bam_size})')

    def _get_bam_file_size(self):
        try:
            size = os.path.getsize(self.bam_file)
            return f"{size / (1024*1024):.1f}MB"
        except:
            return "unknown"

    def bam_analysis_run(self):
        self.printout('subtitle', 'BAM Analysis')

        try:
            import subprocess

            stats_cmd = ['samtools', 'flagstat', self.bam_file]

            if hasattr(self, 'transrate_logger') and self.transrate_logger:
                from core.utils.logging import LoggingSubprocess
                logging_subprocess         = LoggingSubprocess(self.transrate_logger, 'samtools_flagstat')
                returncode, stdout, stderr = logging_subprocess.run_with_logging(stats_cmd, "flagstat")
            else:
                result     = subprocess.run(stats_cmd, capture_output=True, text=True)
                returncode = result.returncode
                stdout     = result.stdout
                stderr     = result.stderr

            if returncode == 0:
                self._parse_bam_stats(stdout)
            else:
                self.printout('warning', f'Failed to analyze BAM file: {stderr}')

        except FileNotFoundError:
            self.printout('warning', 'samtools not found - skipping BAM analysis')
        except Exception as e:
            self.printout('warning', f'Error during BAM analysis: {str(e)}')

    def _parse_bam_stats(self, flagstat_output):
        try:
            lines = flagstat_output.strip().split('\n')
            stats = {}

            for line in lines:
                if 'total' in line and 'secondary' not in line and 'supplementary' not in line:
                    stats['Total Reads'] = int(line.split()[0])
                elif 'mapped' in line and '%' in line:
                    parts = line.split()
                    stats['Mapped Reads'] = int(parts[0])
                    for part in parts:
                        if '%' in part and '(' in line:
                            stats['Mapping Rate (%)'] = float(part.strip('()%'))
                            break
                elif 'properly paired' in line:
                    stats['Properly Paired'] = int(line.split()[0])
                elif 'duplicates' in line:
                    stats['Duplicates'] = int(line.split()[0])

            if stats:
                self.printClass.start_section("BAM Stats")
                self.printClass.add_stage_to_section(stats)
                self.printClass.complete_section()

                if hasattr(self, 'transrate_logger') and self.transrate_logger:
                    self.transrate_logger.log_progress(f"BAM stats extracted: {stats}")

        except Exception as e:
            self.printout('warning', f'Error parsing BAM statistics: {str(e)}')

    def assembly_run(self):
        self.printout('subtitle', 'Assembly Stats')
        assembly_data = {
            'mode'              : self.mode,
            'threads'           : self.threads,
            'assembly'          : self.assembly_file,
            'assembly_base_name': self.assembly_base_name,
            'assemblyDF'        : self.assemblyDF,
            'contigDF'          : self.contigDF,
            'refList'           : self.refList,
            'readCount'         : self.readCount,
            'refCount'          : self.refCount,
            'reference'         : self.reference_file,
            'aHeaders'          : self.aHeaders,
            'dict_file'         : self.dict_file,
            'dict_dir'          : self.dict_dir
        }

        if self.mode > 0:
            assembly_data.update({
                'salmonQuant': os.path.join(self.dict_dir['temp_salmon'], 'quant.sf'),
                'contigCSV'  : os.path.join(self.dict_dir['results']    , f'{self.assembly_base_name}.contigs.csv'),
                'goodContig' : os.path.join(self.dict_dir['results']    , f'good.{self.assembly_base_name}.fa'),
                'badContig'  : os.path.join(self.dict_dir['results']    , f'bad.{self.assembly_base_name}.fa'),
                'scoreOptCSV': os.path.join(self.dict_dir['results']    , 'assembly_score_optimisation.csv')
            })

        assembly_data['assemblyTmp'] = os.path.join(self.dict_dir['temp_analysis'], f'{self.assembly_base_name}.assembly.json')

        AssemblyHub(assembly_data, self.printClass).run(assembly_data)

        self.assemblyDF = assembly_data['assemblyDF']
        self.contigDF   = assembly_data['contigDF']

    def contig_run(self):
        self.printout('subtitle', 'Contig Stats')

        contig_data = {
            'mode'     : self.mode,
            'threads'  : self.threads,
            'contigDF' : self.contigDF,
            'basesDct' : self.basesDct,
            'dict_file': self.dict_file
        }

        ContigHub(contig_data, self.printClass).run(contig_data)

        self.contigDF  = contig_data['contigDF']
        self.basesDct  = contig_data['basesDct']
        self.refList   = contig_data['refList']
        self.refCount  = contig_data['refCount']
        self.readCount = contig_data['readCount']
        print()

        return

    def reference_run(self):
        self.printout('subtitle', 'Reference Analysis')

        from core.assembly.reference import Reference
        multi_mode = getattr(self, 'multi_assembly', False)

        assembly_csv_path = os.path.join(self.dict_dir['results'], 'assembly.csv')

        reference_data = {
            'assembly'          : self.assembly_file,
            'reference'         : self.reference_file,
            'assembly_base_name': self.assembly_base_name,
            'dict_dir'          : self.dict_dir,
            'threads'           : self.threads,
            'multi'             : multi_mode,
            'assemblycsv'       : assembly_csv_path
        }

        self.printClass.start_progress_stage("Reference Analysis")
        self.printClass.update_progress_bar(1, 3, "Setup")

        reference_analyzer = Reference()
        comp_stats         = reference_analyzer.mainRun(reference_data)

        self.printClass.update_progress_bar(2, 3, "Processing results")

        if 'assembly' in comp_stats:
            assembly_stats = comp_stats['assembly']
            for key, value in assembly_stats.items():
                if key in self.assemblyDF.columns:
                    self.assemblyDF.loc[0, key] = value

        self.printClass.update_progress_bar(3, 3, "Complete")
        self.printClass.complete_progress_stage()

        if 'assembly' in comp_stats:
            self.printClass.start_section("Reference Stats")
            self.printClass.add_stage_to_section(comp_stats['assembly'])
            self.printClass.complete_section()

        return

    def get_assembly_results(self):
        results = {}
        if hasattr(self, 'assemblyDF') and self.assemblyDF is not None:
            results['assembly'] = self.assemblyDF.to_dict('records')[0] if len(self.assemblyDF) > 0 else {}
        if hasattr(self, 'contigDF') and self.contigDF is not None:
            results['contigs'] = self.contigDF.to_dict('records')
        return results

    def set_multi_assembly_mode(self, multi_mode: bool = True):
        self.multi_assembly = multi_mode

    def process_csv_files(self):
        assembly_data = {
            'assembly'  : self.assembly_file,
            'assemblyDF': self.assemblyDF,
            'contigDF'  : self.contigDF,
            'cHeaders'  : self.cHeaders,
            'aHeaders'  : self.aHeaders,
            'dict_dir'  : self.dict_dir,
            'contigCSV' : os.path.join(self.dict_dir['results'], f'{self.assembly_base_name}.contigs.csv')
        }

        csv_processor = CSV(self.printClass)
        self.printout('subtitle', 'CSV Processing')

        if self.mode > 0 and self.mode < 4:
            csv_processor.contigCSV(assembly_data)
        if self.assembly_file:
            csv_processor.assemblyCSV(assembly_data)

        return

    def _cleanup_temp_files(self):
        import shutil

        self.transrate_logger.log_stage_start("Cleanup", {"clutter_flag": True})

        try:
            temp_dir = self.dict_dir['temp']
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
                self.transrate_logger.log_file_operation("REMOVE_DIRECTORY", str(temp_dir), True)
                self.printout('info', f'Cleaned up temporary files in {temp_dir}')
            else:
                self.transrate_logger.log_progress("No temp directory found to clean up")

        except Exception as e:
            error_msg = f"Failed to cleanup temp files: {str(e)}"
            self.transrate_logger.log_progress(error_msg, 'error')
            self.printout('warning', error_msg)
