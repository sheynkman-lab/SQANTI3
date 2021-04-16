#%%
import pandas as pd 
import numpy as np

classification = pd.read_table('sqanti_protein_classification_info_chr22.tsv')
# classification = pd.read_csv('../../sqanti_prot_cmp_result_chr22_man_annot_v3.csv')
#%%
classification = (
    classification
        .assign(
            known_nterm=lambda x: x.pr_nterm_diff==0,
            known_cterm=lambda x: x.pr_cterm_diff==0,
            known_nterm_gene=lambda x: x.pr_nterm_gene_diff==0,
            known_cterm_gene=lambda x: x.pr_cterm_gene_diff==0,
            )
)

classification_breakdown = (
    classification
        .groupby(['pr_cat','pr_subcat','known_nterm','known_cterm','known_nterm_gene','known_cterm_gene',])
        .size()
        .reset_index()
        .set_index(['pr_cat','pr_subcat','known_nterm','known_cterm','known_nterm_gene','known_cterm_gene',])
        
)
#%%
#=IF(AND(row['pr_cat']='full-splice_match',OR(row['pr_subcat']='multi-exon',row['pr_subcat']='mono-exon'),row['pr_nterm_diff']=0,row['pr_cterm_diff']=0),'pFSM', 
def is_pfsm(row):
    return (
        row['pr_cat']=='full-splice_match' and 
        (row['pr_subcat']=='multi-exon' or row['pr_subcat'] =='mono-exon') and
        row['pr_nterm_diff']==0 and
        row['pr_cterm_diff']==0
        )

# IF(AND(row['pr_cat']='novel_in_catalog',row['pr_nterm_gene_diff']=0,row['pr_cterm_gene_diff']=0),'pNIC_known_nterm_combo_splice_known_cterm',
def is_pnic_known_nterm_combo_splice_known_cterm(row):
    return(
        row['pr_cat']=='novel_in_catalog' and
        row['pr_nterm_gene_diff']==0 and
        row['pr_cterm_gene_diff']==0
    )

# IF(AND(row['pr_cat']='novel_not_in_catalog',row['pr_nterm_gene_diff']=0,row['pr_cterm_gene_diff']=0),'pNNC_splicing',
def is_pnnc_splicing(row):
    return(
        row['pr_cat']=='novel_not_in_catalog' and
        row['pr_nterm_gene_diff']==0 and
        row['pr_cterm_gene_diff']==0
    ) 

# IF(AND(row['tx_cat']='incomplete-splice_match',row['pr_cat']='incomplete-splice_match',row['pr_cterm_gene_diff']=0,row['pr_nterm_gene_diff']!=0),'Candidate 3frag/Cfrag',
def is_candidate_3frag_cfrag(row):
    return(
        row['tx_cat']=='incomplete-splice_match' and 
        row['pr_cat']=='incomplete-splice_match' and 
        row['pr_nterm_gene_diff']==0 and
        row['pr_cterm_gene_diff'] != 0
    )

# IF(AND(row['pr_cat']='full-splice_match',row['pr_subcat']='multi-exon',OR(row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0)),'pFSM_w_issue',
def is_pfsm_with_issue(row):
    return (
        row['pr_cat']=='full-splice_match' and 
        row['pr_subcat']=='multi-exon' and 
        (row['pr_nterm_diff'] != 0 or row['pr_cterm_diff'] != 0 )
        )
# IF(AND(B2='novel_in_catalog',row['pr_cat']='incomplete-splice_match',D2='combination_of_known_splicesites',row['pr_subcat']='3prime_fragment'),'Legit_new_Nterm?',
def is_possible_new_nterm(row):
    return(
        row['tx_cat']=='novel_in_catalog' and 
        row['pr_cat']=='incomplete-splice_match' and 
        row['tx_subcat']=='combination_of_known_splicesites' and 
        row['pr_subcat']== '3prime_fragment'
    )
        

# IF(AND(OR(row['pr_cat']='intergenic',row['pr_cat']='antisense'),row['pr_subcat']='mono-exon'),'UTR_monoexon_transcript',
def is_utr_monoexon_transcript(row):
    return (
        (row['pr_cat']=='intergenic' or row['pr_cat']=='antisense') and 
        row['pr_subcat']=='mono-exon'
    )

# IF(AND(row['pr_cat']='incomplete-splice_match',row['pr_nterm_gene_diff']=0,row['pr_cterm_gene_diff']!=0),'trunc/altcterm/nmd',
def is_trunc_altcterm_nmd(row):
    return(
        row['pr_cat']=='incomplete-splice_match' and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']!=0

    )
# IF(AND(row['pr_cat']='novel_in_catalog',OR(row['pr_subcat']='mono-exon',row['pr_subcat']='mono-exon_by_intron_retention'),row['pr_nterm_gene_diff']=0),'first_exon_orf',
def is_first_exon_orf(row):
    return(
        row['pr_cat']=='novel_in_catalog' and 
        (row['pr_subcat']=='mono-exon' or row['pr_subcat']=='mono-exon_by_intron_retention') and 
        row['pr_nterm_gene_diff']==0
    )
# IF(AND(row['tx_cat']='novel_in_catalog',row['pr_cat']='novel_in_catalog',OR(row['pr_subcat']='combination_of_known_splicesites',row['pr_subcat']='combination_of_known_junctions',row['pr_subcat']='intron_retention'),row['pr_nterm_gene_diff']=0,row['pr_cterm_gene_diff']!=0),'pNIC_novel_combo_of_known_splice_and_novel_cterm',
def is_pnic_novel_combo_of_known_splice_and_novel_cterm(row):
    return(
        row['tx_cat']=='novel_in_catalog' and 
        row['pr_cat']=='novel_in_catalog' and 
        (
            row['pr_subcat']=='combination_of_known_splicesites' or 
            row['pr_subcat']=='combination_of_known_junctions' or 
            row['pr_subcat']=='intron_retention'
        ) and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']!=0
        
    )
# IF(AND(row['pr_cat']=='novel_not_in_catalog',row['pr_nterm_gene_diff']==0,row['pr_cterm_gene_diff']!=0),'pNNC_novel_splicing_and_cterm',
def is_pnnc_novel_splicing_and_cterm(row):
    return(
        row['pr_cat']=='novel_not_in_catalog' and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']!=0
        
    )
# IF(AND(row['tx_cat']=='incomplete-splice_match',row['pr_cat']=='incomplete-splice_match',row['pr_subcat']=='mono-exon',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']==0),'candidate_3frag_cfrag_monoexon',
def is_candidate_3frag_cfrag_monoexon(row):
    return(
        row['tx_cat']=='incomplete-splice_match' and 
        row['pr_cat']=='incomplete-splice_match' and 
        row['pr_subcat']=='mono-exon' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']==0
    )
# IF(AND(row['pr_cat']=='novel_not_in_catalog',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']==0),'pNNC_novel_nterm_and_splicing',
def is_pNNC_novel_nterm_and_splicing(row):
    return(
        row['pr_cat']=='novel_not_in_catalog' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']==0
    )

# IF(AND(row['pr_cat']=='incomplete-splice_match',row['pr_subcat']=='internal_fragment',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0),'candidate_internal_fragment',
def is_candidate_internal_fragment(row):
    return(
        row['pr_cat']=='incomplete-splice_match' and
        row['pr_subcat']=='internal_fragment' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']!=0
        
    )
# IF(AND(row['pr_subcat']=='mono-exon',row['pr_nterm_gene_diff']==0,row['pr_cterm_gene_diff']==0),'monoexon_w_known_nterm_and_cterm',
def is_monoexon_w_known_nterm_and_cterm(row):
    return(
        row['pr_subcat']=='mono-exon' and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']==0
        
    )
# IF(AND(row['pr_subcat']=='mono-exon',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0),'monoexon_w_unknown_ends',
def is_monoexon_w_unknown_ends(row):
    return(
        row['pr_subcat']=='mono-exon' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']!=0
    )
# IF(AND(row['pr_subcat']=='mono-exon',row['pr_nterm_gene_diff']==0,row['pr_cterm_gene_diff']!=0),'trunc/altcterm/nmd_monoexon',
def is_trunc_altcterm_nmd_monoexon(row):
    return(
        row['pr_subcat']=='mono-exon' and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']!=0
    )
# IF(AND(row['pr_cat']=='incomplete-splice_match',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0),'trunc/novelends',
def is_trunc_novelends(row):
    return(
        row['pr_cat']=='incomplete-splice_match' and
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']!=0
    )
# IF(AND(row['pr_cat']=='incomplete-splice_match',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']==0),'trunc/altnterm',
def is_trunc_altnterm(row):
    return(
        row['pr_cat']=='incomplete-splice_match' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']==0
    )
# IF(AND(row['pr_cat']=='novel_in_catalog',OR(row['pr_subcat']=='combination_of_known_splicesites',row['pr_subcat']=='combination_of_known_junctions'),row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']==0),'pNNC_novel_nterm_known_splice_novel_cterm',
def is_pNNC_novel_nterm_known_splice_novel_cterm(row):
    return(
        row['pr_cat']=='novel_in_catalog' and 
        (
            row['pr_subcat']=='combination_of_known_splicesites' or 
            row['pr_subcat']=='combination_of_known_junctions'
        ) and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']==0
    )
# IF(AND(OR(row['pr_cat']=='full-splice_match',row['pr_cat']=='novel_in_catalog'),row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']==0),'trunc/altnterm_monoexon',
def is_trunc_altnterm_monoexon(row):
    return(
        (
            row['pr_cat']=='full-splice_match' or 
            row['pr_cat']=='novel_in_catalog'
        ) and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']==0
        
    )
# IF(AND(row['pr_cat']=='novel_not_in_catalog',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0),'pNNC_novel_nterm_novel_splice_novel_cterm',
def is_pNNC_novel_nterm_novel_splice_novel_cterm(row):
    return(
        row['pr_cat']=='novel_not_in_catalog' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']!=0
    )
# IF(AND(row['pr_cat']=='novel_in_catalog',row['pr_nterm_gene_diff']!=0,row['pr_cterm_gene_diff']!=0),'pNNC_novel_nterm_combo_splice_novel_cterm',
def is_pNNC_novel_nterm_combo_splice_novel_cterm(row):
    return(
        row['pr_cat']=='novel_in_catalog' and 
        row['pr_nterm_gene_diff']!=0 and 
        row['pr_cterm_gene_diff']!=0
    )
# IF(row['pr_cat']=='fusion','fusion',
def is_fusion(row):
    return(
        row['pr_cat']=='fusion' 
    )
# IF(row['pr_cat']=='intergenic','intergenic',
def is_intergenic(row):
    return(
       row['pr_cat']=='intergenic' 
    )
# IF(AND(OR(row['pr_cat']=='full-splice_match',row['pr_cat']=='incomplete-splice_match'),OR(row['pr_nterm_diff']!=0,row['pr_cterm_diff']!=0),row['pr_nterm_gene_diff']==0,row['pr_cterm_gene_diff']==0),'pNIC_novel_combo_of_known_nterm_and_cterm',''))))))))))))))))))))))))))
def is_pNIC_novel_combo_of_known_nterm_and_cterm(row):
    return(
        (
            row['pr_cat']=='full-splice_match' or 
            row['pr_cat']=='incomplete-splice_match'
        ) and 
        (
            row['pr_nterm_diff']!=0 or row['pr_cterm_diff']!=0
            ) and 
        row['pr_nterm_gene_diff']==0 and 
        row['pr_cterm_gene_diff']==0
    )



# %%
def set_protein_category(row):
    if is_pfsm(row):
        return 'pFSM' 
    if is_pnic_known_nterm_combo_splice_known_cterm(row):   
        return 'pNIC_known_nterm_combo_splice_known_cterm'
    if is_pnnc_splicing(row):
        return 'pNNC_splicing'
    if is_candidate_3frag_cfrag(row):
        return 'candidate_3frac_cfrag'
    if is_pfsm_with_issue(row):
        return 'pfsm_with_issue'
    if is_possible_new_nterm(row):
        return 'possible_new_nterm'
    if is_utr_monoexon_transcript(row):
        return 'utr_monoexon_transcript'
    if is_trunc_altcterm_nmd(row):
        return 'trunc_altcterm_nmd'
    if is_first_exon_orf(row):
        return 'firt_exon_orf'
    if is_pnic_novel_combo_of_known_splice_and_novel_cterm(row):
        return 'pNIC_novel_combo_of_known_splice_and_novel_cterm'
    if is_pnnc_novel_splicing_and_cterm(row):
        return 'pNNC_novel_splicing_and_cterm'
    if is_candidate_3frag_cfrag_monoexon(row):
        return 'candidate_3frag_cfrag_monoexon'
    if is_pNNC_novel_nterm_and_splicing(row):
        return 'pNNC_novel_nterm_and_splicing'
    if is_candidate_internal_fragment(row):
        return 'candidate_internal_fragment'
    if is_monoexon_w_known_nterm_and_cterm(row):
        return 'monoexon_w_known_nterm_and_cterm'
    if is_monoexon_w_unknown_ends(row):
        return 'monoexon_w_unknown_ends'
    if is_trunc_altcterm_nmd_monoexon(row):
        return 'trunc_altcterm_nmd_monoexon'
    if is_trunc_novelends(row):
        return 'trunc_novelends'
    if is_trunc_altnterm(row):
        return 'trunc_altnterm'
    if is_pNNC_novel_nterm_known_splice_novel_cterm(row):
        return 'pNNC_novel_nterm_known_splice_novel_cterm'
    if is_trunc_altnterm_monoexon(row):
        return 'trunc_altnterm_monoexon'
    if is_pNNC_novel_nterm_novel_splice_novel_cterm(row):
        return 'pNNC_novel_nterm_novel_splice_novel_cterm'
    if is_pNNC_novel_nterm_combo_splice_novel_cterm(row):
        return 'pNNC_novel_nterm_combo_splice_novel_cterm'
    if is_fusion(row):
        return 'fusion'
    if is_intergenic(row):
        return 'intergenic'
    if is_pNIC_novel_combo_of_known_nterm_and_cterm(row):
        return 'pNIC_novel_combo_of_known_nterm_and_cterm'
    return 'NO PROTEIN CATEGORY'
    
    
    
    
    


classification['updated_protein_category'] = classification.apply(set_protein_category, axis = 1)
# %%
classification.to_csv('../../updated_classification.tsv', sep='\t')
# %%
