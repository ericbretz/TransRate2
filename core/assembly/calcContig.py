class CalcContig:

    def __init__(self):
        return
    
    def mainRun(self, args):
        main = args
        contDct = {'assembly': {
                            'pFragmentsMapped'  : 0,
                            'pBasesUncovered'   : 0,
                            'pContigsUncovered' : 0,
                            'pContigsLowcovered': 0,
                            'pContigsUncovbase' : 0,
                            'pSoftclipped'      : 0},
                    'contigs': {ref: {'pBasesCovered': 0,} for ref in main['contigDF']['name']}
        }
        if main['mode'] == 2:
            contDct['assembly'].update({
                            'pGoodMappings'     : 0,
                            'pBadMappings'      : 0,
                            'pContigsSegmented' : 0
                            })
            
        contig_names    = main['contigDF']['name']
        bases_uncovered = main['contigDF']['basesUncovered']
        lengths         = main['contigDF']['length']
        p_bases_covered = 1 - (bases_uncovered / lengths)
        contDct['contigs'].update({ref: {'pBasesCovered': p_bases_covered[i]} for i, ref in enumerate(contig_names)})
        
        if main['mode'] == 2:
            both_mapped = main['contigDF']['bothMapped']
            for i, ref in enumerate(contig_names):
                contDct['contigs'][ref]['bothMapped'] = both_mapped[i]

        assemblyDF = main['assemblyDF']
        read_count = main['readCount']
        ref_count  = main['refCount']
        
        contDct['assembly']['pFragmentsMapped']   = assemblyDF['fragmentsMapped'].iloc[0] / read_count
        contDct['assembly']['pBasesUncovered']    = assemblyDF['basesUncovered'].iloc[0] / assemblyDF['bases'].iloc[0]
        contDct['assembly']['pContigsUncovered']  = assemblyDF['contigsUncovered'].iloc[0] / ref_count
        contDct['assembly']['pContigsLowcovered'] = assemblyDF['contigsLowcovered'].iloc[0] / ref_count
        contDct['assembly']['pContigsUncovbase']  = assemblyDF['contigsUncovBase'].iloc[0] / ref_count
        contDct['assembly']['pSoftclipped']       = assemblyDF['softclipped'].iloc[0] / read_count
        
        if main['mode'] == 2:
            contDct['assembly']['pGoodMappings']     = assemblyDF['goodMappings'].iloc[0] / read_count
            contDct['assembly']['pBadMappings']      = assemblyDF['badMappings'].iloc[0] / read_count
            contDct['assembly']['pContigsSegmented'] = assemblyDF['contigsSegmented'].iloc[0] / ref_count
            contDct['assembly']['pGood']             = assemblyDF['goodMappings'].iloc[0] / read_count
        
        return contDct


